const plot=document.getElementById('plot');
const status=document.getElementById('status');
const ready=()=>!rendering&&Boolean(plot._fullLayout?.xaxis);
const frameHeight=()=>send('streamlit:setFrameHeight',{height:document.documentElement.scrollHeight});
const preview=document.createElement('div');preview.id='preview';
const candidate=document.createElement('div');candidate.id='candidate';
const menu=document.createElement('div');menu.id='menu';
menu.innerHTML='<div id="time"></div><button data-label="PASS">PASS</button><button data-label="FAIL">FAIL</button><button data-label="CANCEL">取消</button>';
let args={},drag=null,pending=null,seq=0,rendering=false,clickTimer=null,peakPress=null,lastPeakClick=null;
const doubleClickMs=550;
const send=(type,extra={})=>parent.postMessage({isStreamlitMessage:true,type,...extra},'*');
function currentViewport(){
 const axes={};
 Object.entries(plot._fullLayout||{}).forEach(([name,axis])=>{
  if(/^[xy]axis\d*$/.test(name)&&Array.isArray(axis.range)){
   axes[name]={range:[...axis.range],autorange:axis.autorange,uirevision:axis.uirevision};
  }
 });
 return {revision:plot.layout?.uirevision,axes};
}
function restoreViewport(layout,viewport){
 if(!viewport||viewport.revision!==layout.uirevision)return;
 Object.entries(viewport.axes||{}).forEach(([name,saved])=>{
  const axis=layout[name];
  if(axis&&axis.uirevision===saved.uirevision){
   axis.range=[...saved.range];axis.autorange=saved.autorange;
  }
 });
}
const emit=value=>send('streamlit:setComponentValue',{value:{...value,viewport:currentViewport(),id:Date.now()+':'+(++seq)},dataType:'json'});
function visibleWaveform(layout){
 const wave=args.waveform;
 if(args.signal!=='ecg'||!wave?.length)return null;
 const range=layout.xaxis?.range||layout.xaxis2?.range;
 if(!range)return null;
 const left=Math.max(0,Math.floor(Math.min(...range)*args.sfreq)-1);
 const right=Math.min(wave.length,Math.ceil(Math.max(...range)*args.sfreq)+2);
 const budget=Math.max(1000,Math.ceil(plot.clientWidth*4));
 const stride=Math.max(1,Math.ceil((right-left)/budget)),indices=[];
 for(let start=left;start<right;start+=stride){
  let low=start,high=start;
  for(let i=start+1;i<Math.min(right,start+stride);i++){
   if(wave[i]<wave[low])low=i;
   if(wave[i]>wave[high])high=i;
  }
  indices.push(Math.min(low,high));
  if(low!==high)indices.push(Math.max(low,high));
 }
 // Keep the exact samples used by the editable R-peak markers.
 for(const index of args.peaks)if(index>=left&&index<right)indices.push(index);
 const ordered=[...new Set(indices)].sort((a,b)=>a-b);
 return {x:ordered.map(i=>i/args.sfreq),y:ordered.map(i=>wave[i])};
}
let waveFrame=null;
function refreshVisibleWaveform(update){
 if(rendering||!Object.keys(update).some(k=>/^xaxis\d*\.(range|autorange)/.test(k)))return;
 cancelAnimationFrame(waveFrame);
 waveFrame=requestAnimationFrame(()=>{
  if(rendering)return;
  const axes=currentViewport().axes,wave=visibleWaveform(axes);
  const index=plot.data.findIndex(trace=>trace.name==='Filtered ECG');
  if(wave&&index>=0)Plotly.restyle(plot,{x:[wave.x],y:[wave.y]},[index]);
 });
}
function editPeak(action,index){
 const peaks=action==='DELETE'?args.peaks.filter(p=>p!==index):[...new Set([...args.peaks,index])].sort((a,b)=>a-b);
 args.peaks=peaks;
 const times=peaks.map(p=>p/args.sfreq),rr=times.slice(1).map((t,i)=>t-times[i]);
 const hr=rr.map(t=>t>0?60/t:null),bad=[];
 const [low,high]=args.hr_limits||[40,180],limit=args.rr_deviation??0.2;
 rr.forEach((value,i)=>{
  const local=rr.slice(Math.max(0,i-2),i+3).sort((a,b)=>a-b),mid=Math.floor(local.length/2);
  const median=local.length%2?local[mid]:(local[mid-1]+local[mid])/2;
  if(hr[i]===null||hr[i]<low||hr[i]>high||Math.abs(value-median)/Math.max(median,Number.EPSILON)>limit)bad.push(i);
 });
 const updates={
  'R peaks':{x:times,y:peaks.map(p=>args.waveform[p])},
  'Heart rate':{x:times.slice(1),y:hr},
  'Abnormal RR':{x:bad.map(i=>times[i+1]),y:bad.map(i=>hr[i])}
 };
 // Update all dependent traces together, before the server reruns the page.
 const traces=plot.data.map(trace=>updates[trace.name]?{...trace,...updates[trace.name]}:trace);
 const layout={...plot.layout};
 restoreViewport(layout,currentViewport());
 Plotly.react(plot,traces,layout,plot._context);
 emit({kind:'peak',action,index});
}
function cancel(){drag=null;pending=null;peakPress=null;lastPeakClick=null;clearTimeout(clickTimer);preview.style.display='none';menu.style.display='none';candidate.style.display='none';plot.style.cursor='';}
function coords(e){
 const l=plot._fullLayout,a=l.xaxis,r=plot.getBoundingClientRect();
 return {l,a,r,x:Math.max(a._offset,Math.min(a._offset+a._length,e.clientX-r.left)),y:e.clientY-r.top};
}
function sampleAt(e){
 const {a,x}=coords(e),n=args.waveform.length;
 return Math.max(0,Math.min(n-1,Math.round(a.p2d(x-a._offset)*args.sfreq)));
}
function markerAt(e){
 const {a,r}=coords(e),yaxis=plot._fullLayout.yaxis;
 const px=e.clientX-r.left,py=e.clientY-r.top;
 if(px<a._offset||px>a._offset+a._length||py<yaxis._offset||py>yaxis._offset+yaxis._length)return null;
 let best=null,dist=Infinity;
 for(const index of args.peaks){
  const x=a._offset+a.d2p(index/args.sfreq), y=yaxis._offset+yaxis.d2p(args.waveform[index]);
  if(x<a._offset||x>a._offset+a._length||y<yaxis._offset||y>yaxis._offset+yaxis._length)continue;
  const d=Math.hypot((px-x)/22,(py-y)/32);
  if(d<dist){best=index;dist=d;}
 }
 return dist<=1?best:null;
}
menu.addEventListener('pointerdown',e=>e.stopPropagation());
menu.addEventListener('click',e=>{
 const label=e.target.dataset.label;if(!label)return;
 const bounds=pending;cancel();
 if(label!=='CANCEL'&&bounds)emit({kind:'interval',bounds,label});
});
plot.addEventListener('pointerdown',e=>{
 if(!ready()||e.button!==0||e.target.closest('.modebar')||e.target.closest('#menu'))return;
 if(e.shiftKey)return;
 const {l,a,r,x,y}=coords(e);
 if(e.clientX-r.left<a._offset||e.clientX-r.left>a._offset+a._length||y<l._size.t)return;
 if(args.mode==='peak'){
  if(y<l.yaxis._offset||y>l.yaxis._offset+l.yaxis._length)return;
  peakPress={id:e.pointerId,x:e.clientX,y:e.clientY};
  e.preventDefault();e.stopImmediatePropagation();return;
 }
 cancel();e.preventDefault();e.stopImmediatePropagation();
 drag={id:e.pointerId,x,a,r,label:null};
 plot.setPointerCapture(e.pointerId);
 preview.style.top=l._size.t+'px';preview.style.height=l._size.h+'px';
},{capture:true});
plot.addEventListener('pointermove',e=>{
 if(!ready())return;
 if(drag){
  e.preventDefault();e.stopImmediatePropagation();
  const {x}=coords(e);preview.style.display='block';
  preview.style.left=Math.min(drag.x,x)+'px';preview.style.width=Math.abs(x-drag.x)+'px';return;
 }
 if(args.mode==='peak'&&!e.target.closest('.modebar')){
  const {a,l,x,y}=coords(e),ya=l.yaxis;
  if(y<ya._offset||y>ya._offset+ya._length){candidate.style.display='none';return;}
  const existing=markerAt(e),index=existing===null?sampleAt(e):existing;
  candidate.style.width=candidate.style.height=existing===null?'9px':'18px';
  candidate.style.background=existing===null?'#ef4444':'rgba(239,68,68,.2)';
  candidate.style.border=existing===null?'1px solid white':'2px solid #dc2626';
  plot.style.cursor=existing===null?'crosshair':'pointer';
  candidate.style.display='block';
  candidate.style.left=(a._offset+a.d2p(index/args.sfreq))+'px';
  candidate.style.top=(ya._offset+ya.d2p(args.waveform[index]))+'px';
 }
},{capture:true});
plot.addEventListener('pointerup',e=>{
 if(args.mode==='peak'&&peakPress){
  const press=peakPress;peakPress=null;
  if(rendering||e.shiftKey||e.pointerId!==press.id||Math.hypot(e.clientX-press.x,e.clientY-press.y)>6)return;
  e.preventDefault();e.stopImmediatePropagation();
  const existing=markerAt(e),now=performance.now();
  const previous=lastPeakClick;
  clearTimeout(clickTimer);
  if(existing!==null&&previous&&previous.index===existing&&now-previous.time<=doubleClickMs){
   lastPeakClick=null;candidate.style.display='none';
   editPeak('DELETE',existing);return;
  }
  lastPeakClick={index:existing,time:now};
  if(existing===null){
   const index=sampleAt(e);
   clickTimer=setTimeout(()=>{lastPeakClick=null;editPeak('ADD',index);},doubleClickMs);
  }
  return;
 }
 if(!drag)return;
 e.preventDefault();e.stopImmediatePropagation();
 const {x,r}=coords(e),old=drag;drag=null;
 if(plot.hasPointerCapture(e.pointerId))plot.releasePointerCapture(e.pointerId);
 if(Math.abs(old.x-x)<3){cancel();return;}
 pending=[old.a.p2d(old.x-old.a._offset),old.a.p2d(x-old.a._offset)].sort((a,b)=>a-b);
 menu.querySelector('#time').textContent=pending.map(x=>x.toFixed(3)).join(' – ')+' s';
 menu.style.left=Math.max(0,Math.min(r.width-260,x))+'px';
 menu.style.top=Math.max(0,Math.min(570,e.clientY-r.top))+'px';
 menu.style.display='block';
 menu.querySelector('button').focus({preventScroll:true});
},{capture:true});
plot.addEventListener('click',e=>{
 if(args.mode!=='peak'||e.shiftKey||e.target.closest('.modebar')||e.target.closest('#menu'))return;
 e.preventDefault();e.stopImmediatePropagation();
},{capture:true});
plot.addEventListener('dblclick',e=>{
 if(args.mode!=='peak'||e.target.closest('.modebar'))return;
 e.preventDefault();e.stopImmediatePropagation();
},{capture:true});
plot.addEventListener('pointercancel',cancel);
plot.addEventListener('pointerleave',()=>{peakPress=null;candidate.style.display='none';plot.style.cursor='';});
window.addEventListener('blur',()=>{clearTimeout(clickTimer);cancel();});
window.addEventListener('keydown',e=>{if(e.key==='Escape')cancel();});
plot.addEventListener('wheel',e=>{
 if(e.ctrlKey||!ready())return;
 let dx=e.deltaX||(e.shiftKey?e.deltaY:0);
 if(Math.abs(dx)<.5||(!e.shiftKey&&Math.abs(dx)<Math.abs(e.deltaY)*.35))return;
 e.preventDefault();e.stopImmediatePropagation();cancel();
 const a=plot._fullLayout.xaxis,factor=e.deltaMode===1?16:e.deltaMode===2?a._length:1;
 const shift=dx*factor*(a.range[1]-a.range[0])/a._length,update={};
 Object.keys(plot._fullLayout).filter(k=>/^xaxis\d*$/.test(k)).forEach(k=>{
  const r=plot._fullLayout[k].range;update[k+'.range']=[r[0]+shift,r[1]+shift];
 });Plotly.relayout(plot,update);
},{capture:true,passive:false});
let nextRender=null, renderLoop=false;
async function renderLatest(){
 if(renderLoop)return;
 renderLoop=true;
 try{
  while(nextRender){
   const nextArgs=nextRender;nextRender=null;
   try{
 const viewport=plot._fullLayout?currentViewport():nextArgs.viewport;
 args=nextArgs;clearTimeout(clickTimer);cancel();rendering=true;
 const fig=JSON.parse(args.figure);fig.layout.clickmode='event';fig.layout.dragmode='pan';
 restoreViewport(fig.layout,viewport);
 // Matched axes must start with the same range, including on the first render.
 if(!fig.layout.xaxis.range&&fig.layout.xaxis2?.range)fig.layout.xaxis.range=[...fig.layout.xaxis2.range];
 const wave=visibleWaveform(fig.layout);
 if(wave){const trace=fig.data.find(trace=>trace.name==='Filtered ECG');if(trace)Object.assign(trace,wave);}
 await Plotly.react(plot,fig.data,fig.layout,{scrollZoom:true,displaylogo:false,doubleClick:false,modeBarButtonsToRemove:['select2d','lasso2d']});
 if(!plot._waveViewportListener){plot.on('plotly_relayout',refreshVisibleWaveform);plot._waveViewportListener=true;}
 plot.append(preview,candidate,menu);rendering=false;
    status.hidden=true;
   }catch(error){
    cancel();
    status.textContent='图表加载失败，请刷新页面重试。';status.hidden=false;
    console.error('Chart render failed',error);
   }finally{rendering=false;frameHeight();}
  }
 }finally{renderLoop=false;}
}
window.addEventListener('message',e=>{
 if(e.source!==parent||e.data?.type!=='streamlit:render'||!e.data.args)return;
 nextRender=e.data.args;
 void renderLatest();
});
let resizeFrame=null;
new ResizeObserver(()=>{
 cancelAnimationFrame(resizeFrame);
 resizeFrame=requestAnimationFrame(()=>{
  if(ready())Plotly.Plots.resize(plot);
  frameHeight();
 });
}).observe(document.body);
send('streamlit:componentReady',{apiVersion:1});

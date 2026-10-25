"""Read physical-unit evidence without importing a signal-processing stack."""
from pathlib import Path


def bdf_channel_header(path, channel):
    with Path(path).open('rb') as handle:
        fixed=handle.read(256)
        count=int(fixed[252:256])
        if count <= 0:
            raise ValueError('Invalid BDF channel count: ' + str(path))
        fields={}
        for name,width in [('label',16),('transducer',80),('unit',8),
                           ('physical_min',8),('physical_max',8),
                           ('digital_min',8),('digital_max',8)]:
            block=handle.read(count*width)
            if len(block) != count*width:
                raise ValueError('Truncated BDF header: ' + str(path))
            fields[name]=[block[i*width:(i+1)*width].decode('latin1').strip()
                          for i in range(count)]
    index=fields['label'].index(channel)
    return {key:values[index] for key,values in fields.items()}

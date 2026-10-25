"""Simple interactive launcher: python app.py."""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def prepare(experiment):
    """Stop on failure and leave finalization until after human review."""
    for stage in ('preprocess', 'qc'):
        subprocess.run([sys.executable, str(ROOT/'main.py'), '--experiment', experiment,
                        '--stage', stage], cwd=ROOT, check=True)


def main():
    import project_config as cfg
    from scripts.start_review import review_command

    parser = argparse.ArgumentParser(description='Physio AutoQC launcher')
    parser.add_argument('--experiment', choices=tuple(cfg.EXPERIMENTS))
    parser.add_argument('--action', choices=('prepare', 'review'))
    args = parser.parse_args()
    experiment = args.experiment
    if experiment is None:
        names = list(cfg.EXPERIMENTS)
        print('\n选择实验 / Select experiment:')
        for index, name in enumerate(names, 1):
            print(f'  {index}. {name}')
        while True:
            value = input('输入序号（回车使用默认实验）: ').strip()
            if not value:
                experiment = cfg.EXPERIMENT_NAME
                break
            if value.isdigit() and 1 <= int(value) <= len(names):
                experiment = names[int(value)-1]
                break
            print('请输入列表中的序号。')
    action = args.action
    if action is None:
        print('\n1. 预处理 + 自动质检，然后打开人工复核\n2. 直接打开人工复核（已有处理结果）')
        while True:
            value = input('请选择 [1/2]: ').strip()
            if value in ('1', '2'):
                action = 'prepare' if value == '1' else 'review'
                break
    if action == 'prepare':
        prepare(experiment)
    command, environment, port = review_command(experiment)
    print(f'\n{experiment}: http://127.0.0.1:{port}\n在浏览器打开以上地址；Ctrl+C 停止服务。', flush=True)
    return subprocess.call(command, cwd=ROOT, env=environment)


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except subprocess.CalledProcessError as error:
        print('处理失败，请查看上方错误；未启动后续步骤。', file=sys.stderr)
        raise SystemExit(error.returncode)

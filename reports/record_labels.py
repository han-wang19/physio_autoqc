"""Record labels shared by concert tracks and task-condition experiments."""
import re


def record_parts(record_id):
    match = re.fullmatch(r'A([^_]+)_(.+)', str(record_id))
    if not match:
        raise ValueError('Cannot parse record id: ' + str(record_id))
    condition = match[2]
    music = re.fullmatch(r'music(\d+)', condition)
    return match[1], int(music[1]) if music else condition


def segment_label(value):
    return 'M{}'.format(value) if isinstance(value, int) else str(value).capitalize()


def segment_key(value):
    order = {'resting': 0, 'active': 1, 'resting2': 2, 'passive': 3}
    return (0, value) if isinstance(value, int) else (1, order.get(value, 99), str(value))

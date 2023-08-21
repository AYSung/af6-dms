from pathlib import Path


def set_output_dir(dir_path: str):
    dir_path = Path(dir_path)

    if not dir_path.is_dir():
        dir_path.mkdir()

    return dir_path

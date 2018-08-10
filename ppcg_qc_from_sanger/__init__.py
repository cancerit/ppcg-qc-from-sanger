import subprocess
from typing import Tuple

def exec_subp_and_wait(command: str) -> Tuple[str, int]:
    process = subprocess.Popen(command, shell=True, encoding='UTF-8', stdout=subprocess.PIPE)
    com_return = process.stdout.read()
    returncode = process.wait()
    return com_return.rstrip('\n'), returncode

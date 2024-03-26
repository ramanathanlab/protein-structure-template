from abc import ABC
from enum import Enum


class OutputType(str, Enum):
    stdout = "stdout"
    json = "json"


class Writer(ABC):
    def __init__(self) -> None: ...
    def output(self) -> None: ...

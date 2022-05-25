from .CLI import *
from .src import Normalizer
import os

_ROOT = os.path.abspath(os.path.dirname(__file__))

def get_data(path):
    return os.path.join(_ROOT, 'data', path)

__all__ = [
    'main'
    'Normalizer'
    'get_data'
]
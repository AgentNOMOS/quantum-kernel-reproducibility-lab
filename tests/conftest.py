import os, socket, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def _no_network(*a, **k):
    raise RuntimeError("tests must not use the network")


socket.socket.connect = _no_network
socket.create_connection = _no_network

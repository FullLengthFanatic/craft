"""Smoke test: package imports and exposes __version__."""

import craft


def test_version() -> None:
    assert isinstance(craft.__version__, str)
    assert craft.__version__

from __future__ import annotations

import urllib.request
from collections.abc import Callable
from typing import Any


def fetch_text(
    url: str,
    *,
    timeout: int,
    accept: str = "text/html,text/plain,*/*",
    user_agent: str = "genomi/0.1",
    urlopen: Callable[..., Any] | None = None,
) -> str:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": accept,
            "User-Agent": user_agent,
        },
    )
    open_url = urlopen or urllib.request.urlopen
    with open_url(request, timeout=timeout) as response:
        return response.read().decode("utf-8", "replace")

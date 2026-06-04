from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from typing import Any

JsonObject = dict[str, Any]
_NO_NOT_FOUND_PAYLOAD = object()


def build_api_url(base_url: str, path: str, *, query: dict[str, str] | None = None) -> str:
    encoded_path = "/".join(urllib.parse.quote(part, safe="/") for part in path.split("/"))
    url = base_url.rstrip("/") + "/" + encoded_path.lstrip("/")
    if query:
        url += "?" + urllib.parse.urlencode({key: value for key, value in query.items() if value})
    return url


def fetch_json_with_trace(
    base_url: str,
    path: str,
    *,
    query: dict[str, str] | None = None,
    raw_calls: list[JsonObject],
    timeout: int,
    user_agent: str = "genomi/0.1",
    not_found_payload: Any = _NO_NOT_FOUND_PAYLOAD,
    urlopen: Callable[..., Any] | None = None,
) -> Any:
    """Fetch JSON and append the compact raw-call trace used by API evidence tools."""

    url = build_api_url(base_url, path, query=query)
    call: JsonObject = {"url": url, "status": None, "attempts": 0}
    raw_calls.append(call)
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": user_agent})
    open_url = urlopen or urllib.request.urlopen
    for attempt in range(2):
        call["attempts"] = attempt + 1
        try:
            with open_url(request, timeout=timeout) as response:
                call["status"] = int(getattr(response, "status", 0) or 0)
                call["content_type"] = response.headers.get("content-type")
                body = response.read()
            return json.loads(body.decode("utf-8"))
        except urllib.error.HTTPError as exc:
            call["status"] = exc.code
            if exc.code == 404 and not_found_payload is not _NO_NOT_FOUND_PAYLOAD:
                call["not_found"] = True
                return not_found_payload
            call["error"] = f"HTTP {exc.code}"
            if 400 <= exc.code < 500:
                return None
        except urllib.error.URLError as exc:
            call["error"] = f"URL error: {exc.reason}"
        except TimeoutError as exc:
            call["error"] = f"timeout: {exc}"
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            call["error"] = f"parse error: {exc}"
            return None
        except OSError as exc:
            call["error"] = f"I/O error: {exc}"
        if attempt == 0:
            time.sleep(0.5)
    return None

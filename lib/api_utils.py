#!/usr/bin/env python3
"""
lib/api_utils.py
=================
Shared HTTP client and BioPython utilities used across all modules.

- Throttled, retrying HTTP client for REST APIs
- BioPython-based PDB/sequence helpers
- Common constants (API base URLs, taxonomy, etc.)
"""

import json
import logging
import time
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

log = logging.getLogger(__name__)

# ── API base URLs ────────────────────────────────────────────────
ENSEMBL_REST = "https://rest.ensembl.org"
UNIPROT_REST = "https://rest.uniprot.org"
PDBE_REST = "https://www.ebi.ac.uk/pdbe"
PDBE_GRAPH = "https://www.ebi.ac.uk/pdbe/graph-api"
RCSB_REST = "https://data.rcsb.org"

# ── HTTP helpers ─────────────────────────────────────────────────


def http_get_json(url: str, retries: int = 5, delay: float = 0.5) -> Optional[dict]:
    """GET → parsed JSON with retry & exponential back-off.

    Treats only 404 as a permanent "not found" (returns ``None`` immediately).
    All other HTTP errors and network errors are retried with exponential
    backoff — Ensembl in particular returns transient 400s under load.
    """
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    for attempt in range(retries):
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except HTTPError as e:
            if e.code == 404:
                log.debug("HTTP 404 (not found): %s", url)
                return None
            wait = delay * 2**attempt
            log.warning(
                "HTTP %d on %s – retrying in %.1fs (%d/%d)",
                e.code,
                url,
                wait,
                attempt + 1,
                retries,
            )
            time.sleep(wait)
        except (URLError, TimeoutError, OSError) as e:
            wait = delay * 2**attempt
            log.warning(
                "Network error on %s: %s – retrying in %.1fs (%d/%d)",
                url,
                e,
                wait,
                attempt + 1,
                retries,
            )
            time.sleep(wait)
    log.error("GET failed after %d retries: %s", retries, url)
    return None


def http_post_json(
    url: str,
    payload: dict,
    retries: int = 5,
    delay: float = 0.5,
) -> Optional[dict | list]:
    """POST JSON → parsed JSON with retry & exponential back-off.

    Same retry policy as :func:`http_get_json`: 404 is a permanent miss,
    everything else is retried with exponential backoff.
    """
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    body = json.dumps(payload).encode()
    for attempt in range(retries):
        try:
            req = Request(url, data=body, headers=headers, method="POST")
            with urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode())
        except HTTPError as e:
            if e.code == 404:
                log.debug("HTTP 404 (not found): %s", url)
                return None
            wait = delay * 2**attempt
            log.warning(
                "HTTP %d on %s – retrying in %.1fs (%d/%d)",
                e.code,
                url,
                wait,
                attempt + 1,
                retries,
            )
            time.sleep(wait)
        except (URLError, TimeoutError, OSError) as e:
            wait = delay * 2**attempt
            log.warning(
                "Network error on %s: %s – retrying in %.1fs (%d/%d)",
                url,
                e,
                wait,
                attempt + 1,
                retries,
            )
            time.sleep(wait)
    log.error("POST failed after %d retries: %s", retries, url)
    return None


def http_get_text(
    url: str,
    retries: int = 5,
    delay: float = 0.5,
    accept: str = "*/*",
    timeout: int = 60,
) -> Optional[str]:
    """GET → raw text with retry & exponential back-off.

    Same retry policy as :func:`http_get_json`: 404 is a permanent miss,
    everything else is retried with exponential backoff.
    """
    for attempt in range(retries):
        try:
            req = Request(url, headers={"Accept": accept})
            with urlopen(req, timeout=timeout) as resp:
                return resp.read().decode()
        except HTTPError as e:
            if e.code == 404:
                log.debug("HTTP 404 (not found): %s", url)
                return None
            wait = delay * 2**attempt
            log.warning(
                "HTTP %d on %s – retrying in %.1fs (%d/%d)",
                e.code,
                url,
                wait,
                attempt + 1,
                retries,
            )
            time.sleep(wait)
        except (URLError, TimeoutError, OSError) as e:
            wait = delay * 2**attempt
            log.warning(
                "Network error on %s: %s – retrying in %.1fs (%d/%d)",
                url,
                e,
                wait,
                attempt + 1,
                retries,
            )
            time.sleep(wait)
    log.error("GET failed after %d retries: %s", retries, url)
    return None

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


def http_get_json(url: str, retries: int = 3, delay: float = 0.4) -> Optional[dict]:
    """GET → parsed JSON with retry & rate-limit back-off."""
    headers = {"Accept": "application/json", "Content-Type": "application/json"}
    for attempt in range(retries):
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=60) as resp:
                return json.loads(resp.read().decode())
        except HTTPError as e:
            if e.code == 429:
                wait = delay * 2**attempt
                log.warning("Rate-limited %s – wait %.1fs", url, wait)
                time.sleep(wait)
            elif e.code in (400, 404):
                log.debug("HTTP %d: %s", e.code, url)
                return None
            else:
                log.warning("HTTP %d: %s (attempt %d/%d)", e.code, url, attempt + 1, retries)
                time.sleep(delay)
        except (URLError, TimeoutError, OSError) as e:
            log.warning("Network error %s: %s (%d/%d)", url, e, attempt + 1, retries)
            time.sleep(delay)
    log.error("Failed after %d retries: %s", retries, url)
    return None


def http_get_text(url: str, retries: int = 3, delay: float = 0.4) -> Optional[str]:
    """GET → raw text with retry."""
    for attempt in range(retries):
        try:
            req = Request(url, headers={"Accept": "*/*"})
            with urlopen(req, timeout=60) as resp:
                return resp.read().decode()
        except (HTTPError, URLError, TimeoutError, OSError) as e:
            log.warning("Error %s: %s (%d/%d)", url, e, attempt + 1, retries)
            time.sleep(delay)
    return None

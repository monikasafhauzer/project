"""Self-contained draggable alignment viewer embedded in Streamlit."""
import hashlib
import json
from pathlib import Path


def viewer_html(reference,annotations,rows):
    payload=dict(reference=reference,annotations=annotations,rows=rows)
    fingerprint=hashlib.sha256(json.dumps(payload,sort_keys=True,default=str).encode()).hexdigest()
    payload['fingerprint']=fingerprint
    serialized=json.dumps(payload,default=str).replace('<','\\u003c').replace('>','\\u003e').replace('&','\\u0026')
    return Path(__file__).with_name('viewer.html').read_text().replace('__PAYLOAD__',serialized)

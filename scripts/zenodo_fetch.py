#!/usr/bin/env python3
"""
Download every file from a Zenodo record, verify checksums, and write a manifest.

    python3 zenodo_fetch.py 20584260 -o zenodo_v2
    python3 zenodo_fetch.py 24333082 --manifest-only      # list without downloading

Then, after publishing the new version, diff the manifests to prove nothing is
missing:

    python3 zenodo_fetch.py <v4_id> --manifest-only -o v4
    diff <(cut -d, -f1 zenodo_v2/manifest.csv | sort) \
         <(cut -d, -f1 v4/manifest.csv | sort)

Open records need no token. Standard library only.
"""
import argparse, csv, hashlib, json, os, sys, urllib.request

API = 'https://zenodo.org/api/records/{}'

def get_json(url, token=None):
    req = urllib.request.Request(url, headers={'Accept': 'application/json'})
    if token:
        req.add_header('Authorization', f'Bearer {token}')
    with urllib.request.urlopen(req) as r:
        return json.load(r)

def list_files(rec, record_id):
    """Handle both the legacy Zenodo and the current InvenioRDM response shapes."""
    files = rec.get('files')
    out = []
    if isinstance(files, list):                      # legacy: list of file dicts
        for f in files:
            out.append((f['key'], f.get('size'), f.get('checksum', ''),
                        f['links']['self']))
    elif isinstance(files, dict):                    # InvenioRDM: {"entries": {...}}
        for key, f in (files.get('entries') or {}).items():
            link = (f.get('links') or {}).get('content') \
                or f'https://zenodo.org/api/records/{record_id}/files/{key}/content'
            out.append((key, f.get('size'), f.get('checksum', ''), link))
    return out

def md5(path, buf=1 << 20):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(buf), b''):
            h.update(chunk)
    return h.hexdigest()

def human(n):
    if n is None:
        return '?'
    for u in ('B', 'KB', 'MB', 'GB'):
        if n < 1024:
            return f'{n:.1f}{u}'
        n /= 1024
    return f'{n:.1f}TB'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('record_id', help='Zenodo record id, e.g. 20584260')
    p.add_argument('-o', '--outdir', default=None)
    p.add_argument('--manifest-only', action='store_true')
    p.add_argument('--token', default=os.environ.get('ZENODO_TOKEN'))
    a = p.parse_args()

    outdir = a.outdir or f'zenodo_{a.record_id}'
    os.makedirs(outdir, exist_ok=True)

    rec = get_json(API.format(a.record_id), a.token)
    entries = list_files(rec, a.record_id)
    if not entries:
        sys.exit('No files found. If the record is restricted, set ZENODO_TOKEN.')

    total = sum(s or 0 for _, s, _, _ in entries)
    title = (rec.get('metadata') or {}).get('title', '')
    print(f'Record {a.record_id}: {title[:70]}')
    print(f'{len(entries)} files, {human(total)} total\n')

    rows = []
    for i, (key, size, chk, url) in enumerate(sorted(entries), 1):
        rows.append({'filename': key, 'size_bytes': size, 'checksum': chk})
        if a.manifest_only:
            print(f'  [{i:>3}/{len(entries)}] {key}  ({human(size)})')
            continue
        dest = os.path.join(outdir, key)
        os.makedirs(os.path.dirname(dest) or '.', exist_ok=True)
        if os.path.exists(dest) and size and os.path.getsize(dest) == size:
            print(f'  [{i:>3}/{len(entries)}] {key} — already present, skipping')
        else:
            print(f'  [{i:>3}/{len(entries)}] {key}  ({human(size)}) …', end='', flush=True)
            req = urllib.request.Request(url)
            if a.token:
                req.add_header('Authorization', f'Bearer {a.token}')
            with urllib.request.urlopen(req) as r, open(dest, 'wb') as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            print(' done')
        if chk.startswith('md5:'):
            got = md5(dest)
            if got != chk.split(':', 1)[1]:
                print(f'       !! CHECKSUM MISMATCH for {key}')
            rows[-1]['verified'] = (got == chk.split(':', 1)[1])

    mpath = os.path.join(outdir, 'manifest.csv')
    with open(mpath, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['filename', 'size_bytes', 'checksum', 'verified'])
        w.writeheader()
        for r_ in rows:
            r_.setdefault('verified', '')
            w.writerow(r_)
    print(f'\nManifest: {mpath}')

if __name__ == '__main__':
    main()

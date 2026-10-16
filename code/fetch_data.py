"""Download the velocity cutouts and channel-flow gradients from the Johns
Hopkins Turbulence Databases (JHTDB) with the givernylocal client.

Usage:  JHTDB_TOKEN=<your token> python fetch_data.py --data-dir DATA

A personal JHTDB authorisation token is required (https://turbulence.idies.jhu.edu).
Existing files are skipped, so the script can be restarted after an interruption.
"""
import argparse
import os
import time

import h5py
import numpy as np

import datasets as D


def _client():
    from givernylocal.turbulence_dataset import turb_dataset
    from givernylocal.turbulence_toolkit import getCutout, getData
    return turb_dataset, getCutout, getData


def _retry(fn, attempts=8):
    for k in range(attempts):
        try:
            return fn()
        except Exception:
            if k == attempts - 1:
                raise
            time.sleep(min(1800, 30 * 2 ** k))


def _cutout(ds, getCutout, t, origin):
    x0, y0, z0 = origin
    n = D.CUTOUT
    axes = np.array([[x0, x0 + n - 1], [y0, y0 + n - 1], [z0, z0 + n - 1], [t, t]], dtype=np.int32)
    c = getCutout(ds, "velocity", axes, np.array([1, 1, 1, 1], np.int32), verbose=False)
    da = c[list(c.data_vars)[0]]
    arr = np.asarray(da.transpose("zcoor", "ycoor", "xcoor", "values").values)
    return np.moveaxis(arr.astype(np.float64), -1, 0)          # (3, z, y, x): components u, v, w


def fetch_isotropic(data_dir, token):
    turb_dataset, getCutout, _ = _client()
    for window, w in D.WINDOWS.items():
        ds = None
        for k, (origin, t) in enumerate(zip(w["origins"], w["snapshot"])):
            path = D.cutout_path(data_dir, window, k)
            if os.path.exists(path):
                continue
            if ds is None:
                ds = _retry(lambda: turb_dataset(dataset_title=w["dataset"], output_path=data_dir, auth_token=token))
            u = _retry(lambda: _cutout(ds, getCutout, t, origin))
            with h5py.File(path + ".tmp", "w") as f:
                for name, i in (("u", 0), ("v", 1), ("w", 2)):
                    f.create_dataset(name, data=u[i], compression="gzip", compression_opts=4)
                f.attrs.update({"dataset": w["dataset"], "Re_lambda": w["Re"], "origin": origin, "snapshot": t,
                                "eta": w["eta"], "nu": w["nu"], "dx": 2 * np.pi / w["N"]})
            os.replace(path + ".tmp", path)
            print("fetched", window, k, flush=True)


def fetch_channel(data_dir, token):
    turb_dataset, _, getData = _client()
    ds = _retry(lambda: turb_dataset(dataset_title="channel", output_path=data_dir, auth_token=token))

    def query(t, pts):
        try:
            r = _retry(lambda: getData(ds, "velocity", t, "none", "fd6noint", "gradient", pts, verbose=False))
            return np.asarray(r).reshape(len(pts), 9)
        except Exception as e:
            if "too many storage chunks" in str(e) and len(pts) > 50:
                h = len(pts) // 2
                return np.concatenate([query(t, pts[:h]), query(t, pts[h:])])
            raise

    for stage in D.CHANNEL_STAGES:
        path = D.channel_path(data_dir, stage)
        for b in range(len(D.YPLUS)):
            for k in range(12):
                key = "b%d_k%02d" % (b, k)
                with h5py.File(path, "a") as f:
                    if key in f:
                        continue
                t, pts = D.channel_points(stage, b, k)
                g = query(t, pts)
                with h5py.File(path, "a") as f:
                    d = f.create_dataset(key, data=g)
                    d.attrs.update({"yplus": D.YPLUS[b], "t": t, "points": pts})
                print("fetched channel", stage, key, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", required=True)
    args = ap.parse_args()
    token = os.environ.get("JHTDB_TOKEN")
    if not token:
        raise SystemExit("set the environment variable JHTDB_TOKEN")
    os.makedirs(args.data_dir, exist_ok=True)
    fetch_isotropic(args.data_dir, token)
    fetch_channel(args.data_dir, token)

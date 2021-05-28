from pathlib import Path


def _download(url: str, path: Path):
    try:
        import ipywidgets
        from tqdm.auto import tqdm
    except ImportError:
        from tqdm import tqdm

    from urllib.request import urlopen, Request

    blocksize = 1024 * 8
    blocknum = 0

    try:
        with urlopen(Request(url, headers={"User-agent": "scanpy-user"})) as resp:
            total = resp.info().get("content-length", None)
            with tqdm(
                unit="B",
                unit_scale=True,
                miniters=1,
                unit_divisor=1024,
                total=total if total is None else int(total),
            ) as t, path.open("wb") as f:
                block = resp.read(blocksize)
                while block:
                    f.write(block)
                    blocknum += 1
                    t.update(len(block))
                    block = resp.read(blocksize)

    except (KeyboardInterrupt, Exception):
        # Make sure file doesn’t exist half-downloaded
        if path.is_file():
            path.unlink()
        raise

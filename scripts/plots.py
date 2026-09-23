import csv, os, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.join(os.path.dirname(__file__), "..")
RES = os.path.join(BASE, "results")
OUT = os.path.abspath(os.path.join(BASE, "..", "final_paper", "figs"))
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "lines.linewidth": 1.2,
    "figure.dpi": 300, "axes.grid": True, "grid.alpha": 0.3,
})

STY = {
    "enum":   ("o", "#7b3294", "Enumeration"),
    "sparse": ("s", "#c2a5cf", "Sparse accum."),
    "matzov": ("^", "#2166ac", "MATZOV hybrid"),
    "dense":  ("D", "#1a9850", "Dense FFT (GJ)"),
    "ours":   ("*", "#d6191b", "F2LD (ours)"),
}


def read_kernel():
    rows = []
    with open(os.path.join(RES, "kernel.csv")) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def fig_kernel():
    rows = read_kernel()
    data = {}
    for r in rows:
        mth = r["method"]
        data.setdefault(mth, []).append((float(r["T"]), float(r["time_s"]), float(r["table_bytes"])))
    for k in data:
        data[k].sort()
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.15))
    for mth in ["enum", "sparse", "matzov", "dense", "ours"]:
        if mth not in data:
            continue
        mk, col, lab = STY[mth]
        T = [d[0] for d in data[mth]]
        t = [d[1] * 1e3 for d in data[mth]]
        b = [d[2] / 1024.0 for d in data[mth]]
        ax[0].plot(T, t, mk + "-", color=col, label=lab, markersize=4)
        ax[1].plot(T, b, mk + "-", color=col, label=lab, markersize=4)
    for a in ax:
        a.set_xscale("log", base=2); a.set_yscale("log")
        a.set_xlabel(r"table size $T=p^{k}$")
    ax[0].set_ylabel("median runtime (ms)")
    ax[1].set_ylabel("working memory (KiB)")
    ax[0].legend(loc="upper left", ncol=1, framealpha=0.9)
    fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(OUT, "fig_kernel.eps"), format="eps", bbox_inches="tight")
    plt.close(fig)


def fig_embedded():
    rows = []
    with open(os.path.join(RES, "e2e.csv")) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    ours = [(int(r["h"]), float(r["table_bytes"])) for r in rows if r["method"] == "ours"]
    ours.sort()
    base = {r["method"]: float(r["table_bytes"]) for r in rows if r["method"] != "ours"}
    fig, ax = plt.subplots(figsize=(3.45, 2.35))
    h = [x[0] for x in ours]
    b = [x[1] for x in ours]
    ax.plot(h, b, "*-", color=STY["ours"][2] and "#d6191b", label="F2LD (ours)", markersize=6)
    for mth, col in [("dense", "#1a9850"), ("matzov", "#2166ac")]:
        if mth in base:
            ax.axhline(base[mth], ls="--", color=col, lw=1.0, label=STY[mth][2] + " (no hints)")
    budgets = [("Cortex-M0 (32 KiB)", 32 * 1024, "#888888"),
               ("Cortex-M4 (512 KiB)", 512 * 1024, "#555555"),
               ("Cortex-M7 (2 MiB)", 2 * 1024 * 1024, "#222222")]
    for name, v, col in budgets:
        ax.axhline(v, ls=":", color=col, lw=0.9)
        ax.text(max(h) * 0.62, v * 1.25, name, fontsize=6, color=col)
    ax.set_yscale("log")
    ax.set_xlabel("number of side-channel hints $h$")
    ax.set_ylabel("distinguisher memory (bytes)")
    ax.legend(loc="upper right", framealpha=0.9)
    fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(OUT, "fig_embedded.eps"), format="eps", bbox_inches="tight")
    plt.close(fig)


def fig_success():
    N = []; sp = []; sr = []
    with open(os.path.join(RES, "success.csv")) as f:
        for r in csv.DictReader(f):
            N.append(float(r["N"])); sp.append(float(r["sr_plain"])); sr.append(float(r["sr_prior"]))
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    ax.plot(N, sp, "D-", color="#1a9850", label="without prior", markersize=4)
    ax.plot(N, sr, "*-", color="#d6191b", label="with prior (ours)", markersize=6)
    ax.set_xscale("log")
    ax.set_xlabel("number of dual vectors $N$")
    ax.set_ylabel("empirical success rate")
    ax.set_ylim(-0.03, 1.03)
    ax.legend(loc="lower right", framealpha=0.9)
    fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(OUT, "fig_success.eps"), format="eps", bbox_inches="tight")
    plt.close(fig)


def fig_estimate():
    sec = []
    with open(os.path.join(RES, "estimate.csv")) as f:
        lines = f.read().split("\n")
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("set,h,"):
            start = i + 1; break
    rows = []
    for ln in lines[start:]:
        if not ln.strip():
            continue
        rows.append(ln.split(","))
    series = {}
    for r in rows:
        series.setdefault(r[0], []).append((int(r[1]), float(r[4])))
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    cols = {"Kyber-512": "#d6191b", "Kyber-768": "#2166ac", "Kyber-1024": "#1a9850"}
    for name, col in cols.items():
        if name not in series:
            continue
        s = sorted(series[name])
        ax.plot([x[0] for x in s], [x[1] for x in s], "o-", color=col, label=name, markersize=3)
    for name, v in [("M4 512 KiB", 512 * 1024), ("M7 2 MiB", 2 * 1024 * 1024)]:
        ax.axhline(v, ls=":", color="#444444", lw=0.9)
        ax.text(2, v * 1.4, name, fontsize=6, color="#444444")
    ax.set_yscale("log")
    ax.set_xlabel("pinned coefficients $h$ (side-channel)")
    ax.set_ylabel("online memory (bytes)")
    ax.legend(loc="upper right", framealpha=0.9)
    fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(OUT, "fig_estimate.eps"), format="eps", bbox_inches="tight")
    plt.close(fig)


def fig_extra():
    N = []; sp = []; sr = []
    with open(os.path.join(RES, "success.csv")) as f:
        for r in csv.DictReader(f):
            N.append(float(r["N"])); sp.append(float(r["sr_plain"])); sr.append(float(r["sr_prior"]))
    with open(os.path.join(RES, "estimate.csv")) as f:
        lines = f.read().split("\n")
    start = next(i + 1 for i, ln in enumerate(lines) if ln.startswith("set,h,"))
    series = {}
    for ln in lines[start:]:
        if not ln.strip():
            continue
        r = ln.split(",")
        series.setdefault(r[0], []).append((int(r[1]), float(r[4])))
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.15))
    ax[0].plot(N, sp, "D-", color="#1a9850", label="without prior", markersize=4)
    ax[0].plot(N, sr, "*-", color="#d6191b", label="with prior (ours)", markersize=6)
    ax[0].set_xscale("log"); ax[0].set_xlabel(r"number of dual vectors $N$")
    ax[0].set_ylabel("empirical success rate"); ax[0].set_ylim(-0.03, 1.03)
    ax[0].legend(loc="lower right", framealpha=0.9)
    cols = {"Kyber-512": "#d6191b", "Kyber-768": "#2166ac", "Kyber-1024": "#1a9850"}
    for name, col in cols.items():
        s = sorted(series[name])
        ax[1].plot([x[0] for x in s], [x[1] for x in s], "o-", color=col, label=name, markersize=3)
    for name, v in [("M4 512 KiB", 512 * 1024), ("M7 2 MiB", 2 * 1024 * 1024)]:
        ax[1].axhline(v, ls=":", color="#444444", lw=0.9)
        ax[1].text(2, v * 1.5, name, fontsize=6, color="#444444")
    ax[1].set_yscale("log"); ax[1].set_xlabel(r"pinned coefficients $h$ (side-channel)")
    ax[1].set_ylabel("online memory (bytes)")
    ax[1].legend(loc="upper right", framealpha=0.9)
    fig.tight_layout(pad=0.4)
    fig.savefig(os.path.join(OUT, "fig_extra.eps"), format="eps", bbox_inches="tight")
    plt.close(fig)


def fig_results():
    rows = read_kernel()
    data = {}
    for r in rows:
        data.setdefault(r["method"], []).append((float(r["T"]), float(r["time_s"]), float(r["table_bytes"])))
    for k in data:
        data[k].sort()
    N = []; sp = []; sr = []
    with open(os.path.join(RES, "success.csv")) as f:
        for r in csv.DictReader(f):
            N.append(float(r["N"])); sp.append(float(r["sr_plain"])); sr.append(float(r["sr_prior"]))
    e2e = []
    base = {}
    with open(os.path.join(RES, "e2e.csv")) as f:
        for r in csv.DictReader(f):
            if r["method"] == "ours":
                e2e.append((int(r["h"]), float(r["table_bytes"])))
            else:
                base[r["method"]] = float(r["table_bytes"])
    e2e.sort()

    fig, axg = plt.subplots(2, 2, figsize=(3.5, 3.15))
    ax = [axg[0][0], axg[0][1], axg[1][0], axg[1][1]]
    order = ["enum", "sparse", "matzov", "dense", "ours"]
    for mth in order:
        if mth not in data:
            continue
        mk, col, lab = STY[mth]
        T = [d[0] for d in data[mth]]
        ax[0].plot(T, [d[1] * 1e3 for d in data[mth]], mk + "-", color=col, label=lab, ms=2.4)
        ax[1].plot(T, [d[2] / 1024.0 for d in data[mth]], mk + "-", color=col, label=lab, ms=2.4)
    for a in (ax[0], ax[1]):
        a.set_xscale("log", base=2); a.set_yscale("log"); a.set_xlabel(r"$T=p^{k}$", labelpad=1)
        a.tick_params(labelsize=5.2)
    ax[0].set_ylabel("runtime (ms)", fontsize=6.5); ax[1].set_ylabel("memory (KiB)", fontsize=6.5)
    ax[0].legend(loc="upper left", ncol=1, fontsize=4.3, framealpha=0.9, handlelength=1.1, borderpad=0.25, labelspacing=0.22)
    ax[0].set_title("(a) time", fontsize=7)
    ax[1].set_title("(b) memory", fontsize=7)

    ax[2].plot(N, sp, "D-", color="#1a9850", label="no prior", ms=2.4)
    ax[2].plot(N, sr, "*-", color="#d6191b", label="prior", ms=3.6)
    ax[2].set_xscale("log"); ax[2].set_xlabel(r"$N$", labelpad=1)
    ax[2].set_ylabel("success", fontsize=6.5); ax[2].set_ylim(-0.03, 1.05); ax[2].tick_params(labelsize=5.2)
    ax[2].legend(loc="lower right", fontsize=5.0, framealpha=0.9, handlelength=1.1, borderpad=0.25, labelspacing=0.22)
    ax[2].set_title("(c) recovery", fontsize=7)

    hh = [x[0] for x in e2e]; bb = [x[1] for x in e2e]
    ax[3].plot(hh, bb, "*-", color="#d6191b", label="F2LD", ms=3.6)
    for mth, col in [("dense", "#1a9850"), ("matzov", "#2166ac")]:
        if mth in base:
            ax[3].axhline(base[mth], ls="--", color=col, lw=1.0, label=STY[mth][2].split()[0])
    for name, v in [("M0", 32 * 1024), ("M4", 512 * 1024), ("M7", 2 * 1024 * 1024)]:
        ax[3].axhline(v, ls=":", color="#555555", lw=0.8)
        ax[3].text(max(hh) * 0.80, v * 1.7, name, fontsize=4.6, color="#555555")
    ax[3].set_yscale("log"); ax[3].set_xlabel(r"hints $h$", labelpad=1); ax[3].set_ylabel("bytes", fontsize=6.5)
    ax[3].tick_params(labelsize=5.2)
    ax[3].legend(loc="upper right", fontsize=4.6, framealpha=0.9, handlelength=1.1, borderpad=0.25, labelspacing=0.22)
    ax[3].set_title("(d) hints", fontsize=7)

    fig.tight_layout(pad=0.25, h_pad=0.5, w_pad=0.5)
    fig.savefig(os.path.join(OUT, "fig_results.eps"), format="eps", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    fig_results()
    print("figures written to", OUT)

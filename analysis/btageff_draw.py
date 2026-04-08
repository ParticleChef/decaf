import os
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from hist import loc
from coffea.util import load

plt.style.use(hep.style.CMS)

merged = load("hists/btageff2024.merged")
tagger = "UParTAK4"
outdir = "btageff2024"

flavours = {
    5: "b",
    4: "c",
    0: "lf",
    6: "top",
}

wps = ["loose", "medium", "tight"]

mcs = merged[tagger].keys()

os.makedirs(outdir, exist_ok=True)

for mc in mcs:
    h = merged[tagger][mc]

    for wp in wps:
        for flav, flav_name in flavours.items():
            h_pass = h[{"wp": loc(wp), "btag": loc("pass"), "flavor": loc(flav)}]
            h_fail = h[{"wp": loc(wp), "btag": loc("fail"), "flavor": loc(flav)}]

            num = h_pass.values(flow=False)
            den = h_pass.values(flow=False) + h_fail.values(flow=False)

            eff = np.divide(
                num,
                den,
                out=np.full_like(num, np.nan, dtype=float),
                where=(den != 0),
            )

            pt_edges = h_pass.axes[0].edges
            eta_edges = h_pass.axes[1].edges

            fig, ax = plt.subplots(figsize=(8, 6))

            mesh = ax.pcolormesh(
                eta_edges,
                pt_edges,
                eff,
                shading="auto",
                vmin=0.0,
                vmax=1.0,
            )

            cbar = fig.colorbar(mesh, ax=ax)
            cbar.set_label("Efficiency")

            hep.cms.label(
                llabel="Simulation",
                rlabel="(13.6 TeV)",
                loc=0,
                ax=ax,
            )

            ax.set_xlabel(r"Jet |$\eta$|")
            ax.set_ylabel(r"Jet $p_{T}$ (GeV)")
            ax.set_ylim(30,1000)
            ax.set_yscale("log")

            eta_centers = 0.5 * (eta_edges[:-1] + eta_edges[1:])
            pt_centers = np.sqrt(pt_edges[:-1] * pt_edges[1:])

            norm = plt.Normalize(vmin=0.0, vmax=1.0)
            cmap = plt.get_cmap("viridis")

            for ipt, ptc in enumerate(pt_centers):
                for ieta, etac in enumerate(eta_centers):
                    val = eff[ipt, ieta]
                    if np.isfinite(val):
                        rgba = cmap(norm(val))
                        r, g, b, _ = rgba
                        luminance = 0.299 * r + 0.587 * g + 0.114 * b
                        text_color = "white" if luminance < 0.5 else "black"

                        ax.text(
                            etac,
                            ptc,
                            f"{val:.3f}",
                            ha="center",
                            va="center",
                            fontsize=14,
                            color=text_color,
                            fontweight="bold",
                        )

            fig.tight_layout()
            fig.savefig(f"{outdir}/{mc}_{flav_name}_{wp}.png", dpi=200)
            plt.close(fig)

            print(f"saved: {outdir}/{mc}_{flav_name}_{wp}.png")
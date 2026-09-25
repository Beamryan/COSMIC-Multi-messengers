import pandas as pd, plotly.graph_objects as go, plotly.io as pio
pio.renderers.default = "browser"


df = pd.read_csv("visible_stars50.csv")
x, y, z, A_V = (df[col].to_numpy() for col in ["xHx(kpc)", "yHx(kpc)", "zHx(kpc)", "extinction"])

# colour map  (Turbo → blue→yellow)
cmin, cmax = A_V.min(), A_V.max()
colors = (A_V - cmin) / (cmax - cmin + 1e-9)    # 0…1 for colourscale

# build line segments Sun→star, with None separators
xs, ys, zs, cs = [], [], [], []
for xi, yi, zi, ci in zip(x, y, z, colors):
    xs += [0, xi, None]
    ys += [0, yi, None]
    zs += [0, zi, None]
    cs += [ci, ci, ci]   # constant colour along the segment

fig = go.Figure()

fig.add_trace(go.Scatter3d(
    x=xs, y=ys, z=zs,
    mode="lines",
    line=dict(width=2,
            colorscale="Turbo",
            cmin=0, cmax=1,
            color=cs,
            reversescale=False),
    name="extinction rays"
))

# optional: dots at star positions
fig.add_trace(go.Scatter3d(
    x=x, y=y, z=z,
    mode="markers",
    marker=dict(size=2,
                color=colors, 
                colorscale="Turbo", 
                cmin=0, cmax=1,
                opacity=0.8),
    name="binary"
))

fig.update_layout(
    scene=dict(aspectmode="data",
               xaxis_title="X [kpc]", yaxis_title="Y [kpc]", zaxis_title="Z [kpc]"),
    title="White-dwarf binaries – extinction coloured sight-lines",
    legend=dict(itemsizing="constant")
)
fig.show()
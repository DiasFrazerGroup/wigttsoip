# Vendored from ../../alphagenome_finetuning_rna/figures/figutils.py (cm + set_figure_style
# only - the rest of that file depends on statannotations, which we don't need here).

cm = 1 / 2.54


def set_figure_style():

    import matplotlib as mpl
    import matplotlib.font_manager as fm
    import seaborn as sns

    # Font settings — use Arial if available, otherwise DejaVu Sans
    available = {f.name for f in fm.fontManager.ttflist}
    mpl.rcParams['font.family'] = 'Arial' if 'Arial' in available else 'DejaVu Sans'
    mpl.rcParams["font.size"] = 6

    mpl.rcParams['xtick.labelsize'] = 6
    mpl.rcParams['ytick.labelsize'] = 6

    mpl.rcParams['axes.labelsize'] = 8
    mpl.rcParams['axes.titlesize'] = 8

    mpl.rcParams['legend.fontsize'] = 6
    mpl.rcParams['legend.title_fontsize'] = 8

    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['figure.figsize'] = (2.5, 2.5)
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.right'] = False

    # thin black outlines on bars/patches, matching the thinner axis lines above
    mpl.rcParams['patch.linewidth'] = 0.5
    mpl.rcParams['patch.edgecolor'] = 'black'
    mpl.rcParams['patch.force_edgecolor'] = True

    mpl.rcParams['figure.dpi'] = 150
    mpl.rcParams['savefig.dpi'] = 300

    sns.set_context(rc={'figure.figsize': (2.5, 2.5)})

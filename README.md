# ClusterPrepForDFT
Working on systematically generating clusters to do solid state modelling in a molecular DFT code, e.g. NWChem.

## Installation

ClusterPrepForDFT requires pymatgen, which can be installed with:
`pip install pymatgen`

Windows: If you run into issues with `spglib` requiring microsoft visual studio build tools, you have two options:
- (1) Install instead through conda: `conda install -c conda-forge spglib`
- (2) Install build tools and run the pip install command again. See (https://stackoverflow.com/questions/64261546/python-cant-install-packages) for more details.

Trial text for git push


Quick command for installing visual studio build tools:
Download build tools directly: https://aka.ms/vs/16/release/vs_buildtools.exe
Then in command prompt or powershell:
`vs_buildtools.exe --norestart --passive --downloadThenInstall --includeRecommended --add Microsoft.VisualStudio.Workload.NativeDesktop --add Microsoft.VisualStudio.Workload.VCTools --add Microsoft.VisualStudio.Workload.MSBuildTools`

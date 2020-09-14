
- [`water-scaled-arrow.pdf`](water-scaled-arrow.pdf): Image depicting molecular interactions between atoms being turned off and on during an NCMC translational water move. [`water-scaled.pdf`](water-scaled.pdf) is the same image, but without an arrow indicating water movement. Made in [keynote](https://www.apple.com/keynote/) by Danielle Bergazin.

- [`3system-crystal-water-circled.png`](3system-crystal-water-circled.png) and [`3system-crystal-water-circled.jpg`](3system-crystal-water-circled.jpg): Image of three of the systems used to test the ability of NCMC/MD water hopping to allow the exchange of water. The crystallographic water molecules are circled. Individual system images were generated with UCSF Chimera and the final image was made with Inkscape by Danielle Bergazin.

- [`radius-example-figure.png`](radius-example-figure.png): Image depicting a user-defined radius that covers a particular area of interest. The MUP-1 protein-ligand system is used in this example. The radius (indicated by the black dashed line) defines a sphere around a user-selected atom (represented by a blue star) in the system.

- [`mup-acceptance-ratio-plot.pdf`](mup-acceptance-ratio-plot.pdf): Plot showing the rate of water transfer from bulk to the internal hydration site in MUP-1 versus increase in NCMC steps. The plot was generated with [Matplotlib](https://matplotlib.org/#) by Danielle Bergazin.

- [`graphene-wall-plot-fig-iter.pdf`](graphene-wall-plot-fig-iter.pdf): Figure showing the water box system with dividing graphene sheets and the densities in the two regions stabilizing over time. The plot was generated with [Matplotlib](https://matplotlib.org/#), the image of the system was generated with UCSF Chimera and the final image was made with Inkscape by Danielle Bergazin. Figure is in `.pdf` and `.svg` format. `graphene-wall-plot-fig-iter.svg` is [Inkscape](https://inkscape.org) compatible.

- [`blues-cycle.pdf`](blues-cycle.pdf): Image depicting the workflow of BLUES water hopping proposals. Made in [keynote](https://www.apple.com/keynote/) by Danielle Bergazin.

- [`MUP-1-acceptance-and-force-evaluation-and-wallclock-time-1250-2500-5000-30000-table.tgn`](MUP-1-acceptance-and-force-evaluation-and-wallclock-time-1250-2500-5000-30000-table.tgn): Table listing the average acceptance rate of all BLUES moves, the average number of force evaluations across 10-12 replicates for the buried cavity in the MUP-1 system to become hydrated, and the average wallclock time in hours for BLUES to hydrate MUP-1. Each simulation was run for 10000 BLUES iterations, where each iteration consisted of a single NCMC move (consisting of n NCMC steps) and 1000 MD steps.

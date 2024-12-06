# waves-currentshear-windstress

Codes used to generate graphics used in manuscript to be submitted to _Environmental Fluid Mechanics_:

"Accounting for Waves and Current Shear in Ocean Wind Stress Parameterization"
 
D. G. Ortiz-Suslow, N. J. M. Laxague, M. Curcic, & J.-V. I. Bjorkqvist, 2024

## Manuscript Figure Generation
Navigate to _2024\_EFM\_manuscript_ and run _EFM\_2024\_figure\_gen\_script.m_.

## Core Functions Associated With Wind-Wave-Current Interaction Framework
Navigate to _core\_functions_.

The script _compute\_form\_drag\_effective\_current\_OOI\_script.m_ is the command center, though it won't work out-of-the-box (in its present state, it requires _.mat_ files which contain the OOI and NDBC data chunks).

The functions called by that "command center" script should run without issue.
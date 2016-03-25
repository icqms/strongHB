# strongHB
2Y3L pdb quantum refinement and optimization

1) pdb2hin   created hin file used as coord and info database in all calcs

2) checkres looks for errros in pdb file, makes sure all res agree with templates (eg, adds hydrogens)

3) opt-res   is control program for splitting optimization (or refinemet) into small jobs; first calcs cam-b3lyp/6-31(+)G*/amber

4) reenter   collects outputs from calcs initiated by opt-res and updates hin file

these calcs are relevant to JACS article 10.1021/jacs.5b11012  on strong hydrogen bondin

readhin.for 
 -reads hin file
 -writes hin file
 -pbc 
 -expand copies
 -add hydrogen atoms
 -delete residues
 -copy residues
 
template.for
  -deals with templates, only used by checkres.for
  
checkres.for
  -check for errors in the pdb.
 
reenter.for
 -bring back all the bits from g09.log files - into an updated hin file.

//example to compile
ifort readhin.obj checkres.for -o checkres

//execute
./checkres 2y3l
//
 this will make a 2y3lc.hin file


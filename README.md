# strongHB
2Y3L pdb quantum refinement and optimization
1) pdb2hin   created hin file used as coord and info database in all calcs
2) checkres looks for errros in pdb file, makes sure all res agree with templates (eg, adds hydrogens)
3) opt-res   is control program for splitting optimization (or refinemet) into small jobs
4) reenter   collects outputs from calcs initiated by opt-res and updates hin file

these calcs are relevant to JACS article 10.1021/jacs.5b11012  on strong hydrogen bondin
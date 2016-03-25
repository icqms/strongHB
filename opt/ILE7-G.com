g09 << +++ > ILE7-G.log
%mem=1900MB
%Nprocshared=4
%Chk=/short/d63/reimers/ILE7-G
#P ONIOM(cam-b3lyp:AMBER=hardfirst)=EmbedCharge iop(4/28=5)
   opt(z-matrix,maxcycle=7,nrscale) nosym scf(conver=6) 

  lysozyme partial opt ILE7-G

  1 1 -1 1 -1 1
 N-N--0.416 0 x0001 y0001 z0001    
 C-CT--0.060 0 x0002 y0002 z0002    
 C-C-0.597 0 x0003 y0003 z0003    
 O-O--0.568 0 x0004 y0004 z0004    
 C-CT-0.130 0 x0005 y0005 z0005    
 C-CT--0.043 0 x0006 y0006 z0006    
 C-CT--0.320 0 x0007 y0007 z0007    
 C-CT--0.066 0 x0008 y0008 z0008    
 H-H-0.272 0 x0009 y0009 z0009    
 H-H1-0.087 0 x0010 y0010 z0010    
 H-HC-0.019 0 x0011 y0011 z0011    
 H-HC-0.024 0 x0012 y0012 z0012    
 H-HC-0.024 0 x0013 y0013 z0013    
 H-HC-0.088 0 x0014 y0014 z0014    
 H-HC-0.088 0 x0015 y0015 z0015    
 H-HC-0.088 0 x0016 y0016 z0016    
 H-HC-0.019 0 x0017 y0017 z0017    
 H-HC-0.019 0 x0018 y0018 z0018    
 H-HC-0.019 0 x0019 y0019 z0019    
 N-N3-0.175 0 x0020 y0020 z0020 M  
 C-CT--0.024 0 x0021 y0021 z0021 M  
 C-C-0.597 0 x0022 y0022 z0022 M  
 O-O--0.568 0 x0023 y0023 z0023 M  
 C-CT-0.034 0 x0024 y0024 z0024 M  
 C-CT-0.002 0 x0025 y0025 z0025 M  
 S-S--0.274 0 x0026 y0026 z0026 M  
 C-CT--0.054 0 x0027 y0027 z0027 M  
 H-H-0.227 0 x0028 y0028 z0028 M  
 H-H1-0.088 0 x0029 y0029 z0029 M  
 H-HC-0.024 0 x0030 y0030 z0030 M  
 H-HC-0.024 0 x0031 y0031 z0031 M  
 H-H1-0.044 0 x0032 y0032 z0032 M  
 H-H1-0.044 0 x0033 y0033 z0033 M  
 H-H1-0.068 0 x0034 y0034 z0034 M  
 H-H1-0.068 0 x0035 y0035 z0035 M  
 H-H1-0.068 0 x0036 y0036 z0036 M  
 H-H-0.227 0 x0037 y0037 z0037 M  
 H-H-0.227 0 x0038 y0038 z0038 M  
 N-N--0.416 0 x0039 y0039 z0039 M  
 C-CT--0.087 0 x0040 y0040 z0040 M  
 C-C-0.597 0 x0041 y0041 z0041 M  
 O-O--0.568 0 x0042 y0042 z0042 M  
 C-CT-0.298 0 x0043 y0043 z0043 M  
 C-CT--0.319 0 x0044 y0044 z0044 M  
 C-CT--0.319 0 x0045 y0045 z0045 M  
 H-H-0.272 0 x0046 y0046 z0046 M  
 H-H1-0.097 0 x0047 y0047 z0047 M  
 H-HC--0.029 0 x0048 y0048 z0048 M  
 H-HC-0.079 0 x0049 y0049 z0049 M  
 H-HC-0.079 0 x0050 y0050 z0050 M  
 H-HC-0.079 0 x0051 y0051 z0051 M  
 H-HC-0.079 0 x0052 y0052 z0052 M  
 H-HC-0.079 0 x0053 y0053 z0053 M  
 H-HC-0.079 0 x0054 y0054 z0054 M  
 N-N--0.416 0 x0055 y0055 z0055 M H-H1--0.019 
 C-CT--0.087 0 x0056 y0056 z0056    
 C-C-0.597 0 x0057 y0057 z0057    
 O-O--0.568 0 x0058 y0058 z0058    
 C-CT-0.298 0 x0059 y0059 z0059 M H-H1--0.019 
 C-CT--0.319 0 x0060 y0060 z0060 M  
 C-CT--0.319 0 x0061 y0061 z0061 M  
 H-H-0.272 0 x0062 y0062 z0062 M  
 H-H1-0.097 0 x0063 y0063 z0063    
 H-HC--0.030 0 x0064 y0064 z0064 M  
 H-HC-0.079 0 x0065 y0065 z0065 M  
 H-HC-0.079 0 x0066 y0066 z0066 M  
 H-HC-0.079 0 x0067 y0067 z0067 M  
 H-HC-0.079 0 x0068 y0068 z0068 M  
 H-HC-0.079 0 x0069 y0069 z0069 M  
 H-HC-0.079 0 x0070 y0070 z0070 M  
 N-N--0.416 0 x0071 y0071 z0071    
 C-CT--0.087 0 x0072 y0072 z0072    
 C-C-0.597 0 x0073 y0073 z0073    
 O-O--0.568 0 x0074 y0074 z0074    
 C-CT-0.298 0 x0075 y0075 z0075    
 C-CT--0.319 0 x0076 y0076 z0076    
 C-CT--0.319 0 x0077 y0077 z0077    
 H-H-0.272 0 x0078 y0078 z0078    
 H-H1-0.097 0 x0079 y0079 z0079    
 H-HC--0.030 0 x0080 y0080 z0080    
 H-HC-0.079 0 x0081 y0081 z0081    
 H-HC-0.079 0 x0082 y0082 z0082    
 H-HC-0.079 0 x0083 y0083 z0083    
 H-HC-0.079 0 x0084 y0084 z0084    
 H-HC-0.079 0 x0085 y0085 z0085    
 H-HC-0.079 0 x0086 y0086 z0086    
 N-N--0.416 0 x0087 y0087 z0087    
 C-CT-0.034 0 x0088 y0088 z0088    
 C-C-0.412 0 x0089 y0089 z0089    
 O-O2--0.568 0 x0090 y0090 z0090    
 C-CT--0.182 0 x0091 y0091 z0091    
 O-O2--0.815 0 x0092 y0092 z0092    
 H-H-0.272 0 x0093 y0093 z0093    
 H-H1-0.082 0 x0094 y0094 z0094    
 H-HC-0.060 0 x0095 y0095 z0095    
 H-HC-0.060 0 x0096 y0096 z0096    
 H-HC-0.060 0 x0097 y0097 z0097    
 N-N--0.416 0 x0098 y0098 z0098 M  
 C-CT--0.025 0 x0099 y0099 z0099 M  
 C-C-0.597 0 x0100 y0100 z0100 M  
 O-O--0.568 0 x0101 y0101 z0101 M  
 H-H-0.272 0 x0102 y0102 z0102 M  
 H-H1-0.070 0 x0103 y0103 z0103 M  
 H-H1-0.070 0 x0104 y0104 z0104 M  
 N-N--0.416 0 x0105 y0105 z0105 M  
 C-CT--0.025 0 x0106 y0106 z0106 M  
 C-C-0.597 0 x0107 y0107 z0107 M  
 O-O--0.568 0 x0108 y0108 z0108 M  
 H-H-0.272 0 x0109 y0109 z0109 M  
 H-H1-0.070 0 x0110 y0110 z0110 M  
 H-H1-0.070 0 x0111 y0111 z0111 M  
 N-N--0.416 0 x0112 y0112 z0112 M  
 C-CT--0.087 0 x0113 y0113 z0113 M  
 C-C-0.546 0 x0114 y0114 z0114 M  
 O-O--0.559 0 x0115 y0115 z0115 M  
 C-CT-0.298 0 x0116 y0116 z0116 M  
 C-CT--0.319 0 x0117 y0117 z0117 M  
 C-CT--0.319 0 x0118 y0118 z0118 M  
 H-H-0.272 0 x0119 y0119 z0119 M  
 H-H1-0.097 0 x0120 y0120 z0120 M  
 H-HC--0.030 0 x0121 y0121 z0121 M  
 H-HC-0.079 0 x0122 y0122 z0122 M  
 H-HC-0.079 0 x0123 y0123 z0123 M  
 H-HC-0.079 0 x0124 y0124 z0124 M  
 H-HC-0.079 0 x0125 y0125 z0125 M  
 H-HC-0.079 0 x0126 y0126 z0126 M  
 H-HC-0.079 0 x0127 y0127 z0127 M  
 N-N3-0.175 0 x0128 y0128 z0128 M H-H1--0.047 
 C-CT--0.024 0 x0129 y0129 z0129    
 C-C-0.597 0 x0130 y0130 z0130    
 O-O--0.568 0 x0131 y0131 z0131    
 C-CT-0.034 0 x0132 y0132 z0132 M H-H1--0.047 
 C-CT-0.002 0 x0133 y0133 z0133 M  
 S-S--0.274 0 x0134 y0134 z0134 M  
 C-CT--0.054 0 x0135 y0135 z0135 M  
 H-H-0.227 0 x0136 y0136 z0136 M  
 H-H1-0.088 0 x0137 y0137 z0137    
 H-HC-0.024 0 x0138 y0138 z0138 M  
 H-HC-0.024 0 x0139 y0139 z0139 M  
 H-H1-0.044 0 x0140 y0140 z0140 M  
 H-H1-0.044 0 x0141 y0141 z0141 M  
 H-H1-0.058 0 x0142 y0142 z0142 M  
 H-H1-0.068 0 x0143 y0143 z0143 M  
 H-H1-0.068 0 x0144 y0144 z0144 M  
 H-H-0.227 0 x0145 y0145 z0145 M  
 H-H-0.227 0 x0146 y0146 z0146 M  
 N-N--0.416 0 x0147 y0147 z0147    
 C-CT--0.087 0 x0148 y0148 z0148    
 C-C-0.597 0 x0149 y0149 z0149    
 O-O--0.568 0 x0150 y0150 z0150    
 C-CT-0.298 0 x0151 y0151 z0151    
 C-CT--0.319 0 x0152 y0152 z0152    
 C-CT--0.319 0 x0153 y0153 z0153    
 H-H-0.272 0 x0154 y0154 z0154    
 H-H1-0.097 0 x0155 y0155 z0155    
 H-HC--0.030 0 x0156 y0156 z0156    
 H-HC-0.079 0 x0157 y0157 z0157    
 H-HC-0.079 0 x0158 y0158 z0158    
 H-HC-0.079 0 x0159 y0159 z0159    
 H-HC-0.079 0 x0160 y0160 z0160    
 H-HC-0.079 0 x0161 y0161 z0161    
 H-HC-0.079 0 x0162 y0162 z0162    
 N-N--0.416 0 x0163 y0163 z0163 M  
 C-CT--0.087 0 x0164 y0164 z0164 M  
 C-C-0.597 0 x0165 y0165 z0165 M  
 O-O--0.568 0 x0166 y0166 z0166 M  
 C-CT-0.292 0 x0167 y0167 z0167 M  
 C-CT--0.319 0 x0168 y0168 z0168 M  
 C-CT--0.133 0 x0169 y0169 z0169 M  
 H-H-0.272 0 x0170 y0170 z0170 M  
 H-H1-0.097 0 x0171 y0171 z0171 M  
 H-HC--0.018 0 x0172 y0172 z0172 M  
 H-HC-0.079 0 x0173 y0173 z0173 M  
 H-HC-0.039 0 x0174 y0174 z0174 M  
 H-HC-0.079 0 x0175 y0175 z0175 M  
 H-HC-0.026 0 x0176 y0176 z0176 M  
 H-HC-0.039 0 x0177 y0177 z0177 M  
 H-HC-0.014 0 x0178 y0178 z0178 M  
 N-N--0.416 0 x0179 y0179 z0179 M  
 C-CT--0.087 0 x0180 y0180 z0180 M  
 C-C-0.437 0 x0181 y0181 z0181 M  
 O-O--0.177 0 x0182 y0182 z0182 M  
 C-CT-0.298 0 x0183 y0183 z0183 M  
 C-CT--0.319 0 x0184 y0184 z0184 M  
 C-CT--0.319 0 x0185 y0185 z0185 M  
 H-H-0.271 0 x0186 y0186 z0186 M  
 H-H1-0.097 0 x0187 y0187 z0187 M  
 H-HC--0.030 0 x0188 y0188 z0188 M  
 H-HC-0.079 0 x0189 y0189 z0189 M  
 H-HC-0.079 0 x0190 y0190 z0190 M  
 H-HC-0.079 0 x0191 y0191 z0191 M  
 H-HC-0.079 0 x0192 y0192 z0192 M  
 H-HC-0.079 0 x0193 y0193 z0193 M  
 H-HC-0.079 0 x0194 y0194 z0194 M  
 C-CT--0.025 0 x0195 y0195 z0195 M  
 H-H-0.272 0 x0196 y0196 z0196 M  
 N-N--0.416 0 x0197 y0197 z0197 M  
 H-C-0.597 0 x0198 y0198 z0198 M  
 H-H1-0.070 0 x0199 y0199 z0199 M  
 H-H1-0.070 0 x0200 y0200 z0200 M  
 C-CT--0.025 0 x0201 y0201 z0201 M  
 O-O--0.568 0 x0202 y0202 z0202 M  
 H-N--0.416 0 x0203 y0203 z0203 M  
 C-C-0.597 0 x0204 y0204 z0204 M  
 H-H1-0.070 0 x0205 y0205 z0205 M  
 H-H1-0.070 0 x0206 y0206 z0206 M  
 C-CT--0.025 0 x0207 y0207 z0207 M  
 O-O--0.568 0 x0208 y0208 z0208 M  
 H-N--0.133 0 x0209 y0209 z0209 M  
 C-C-0.581 0 x0210 y0210 z0210 M  
 H-CT-0.048 0 x0211 y0211 z0211 M  
 H-H1-0.023 0 x0212 y0212 z0212 M  
 C-CT--0.011 0 x0213 y0213 z0213 M  
 H-H-0.087 0 x0214 y0214 z0214 M  
 N-N--0.122 0 x0215 y0215 z0215 M  
 H-C-0.090 0 x0216 y0216 z0216 M  
 H-CT-0.023 0 x0217 y0217 z0217 M  
 H-H1-0.016 0 x0218 y0218 z0218 M  
 C-CT--0.025 0 x0219 y0219 z0219    
 H-H-0.272 0 x0220 y0220 z0220    
 N-N--0.416 0 x0221 y0221 z0221    
 H-C-0.597 0 x0222 y0222 z0222    
 H-H1-0.070 0 x0223 y0223 z0223    
 H-H1-0.070 0 x0224 y0224 z0224    
 C-CT--0.025 0 x0225 y0225 z0225 M  
 O-O--0.553 0 x0226 y0226 z0226 M  
 H-N--0.375 0 x0227 y0227 z0227 M  
 C-C-0.597 0 x0228 y0228 z0228 M  
 H-H1-0.070 0 x0229 y0229 z0229 M  
 H-H1-0.070 0 x0230 y0230 z0230 M  
 C-CT--0.013 0 x0231 y0231 z0231 M  
 H-H-0.212 0 x0232 y0232 z0232 M  
 N-N--0.240 0 x0233 y0233 z0233 M  
 H-C-0.183 0 x0234 y0234 z0234 M  
 H-CT-0.015 0 x0235 y0235 z0235 M  
 H-H1-0.014 0 x0236 y0236 z0236 M  
 Variables:
x0009 =    -3.7550
y0009 =     9.3531
z0009 =    -0.4425
x0010 =    -1.2934
y0010 =    10.6566
z0010 =    -1.6307
x0011 =    -3.7992
y0011 =     9.8702
z0011 =    -3.2381
x0012 =    -2.1753
y0012 =     8.0032
z0012 =    -2.2860
x0013 =    -1.1919
y0013 =     8.5016
z0013 =    -3.7046
x0014 =    -2.7028
y0014 =    11.7833
z0014 =    -4.2099
x0015 =    -2.2228
y0015 =    10.3542
z0015 =    -5.1872
x0016 =    -1.0786
y0016 =    11.0481
z0016 =    -3.9884
x0017 =    -3.8227
y0017 =     7.1027
z0017 =    -3.5686
x0018 =    -2.4278
y0018 =     6.8666
z0018 =    -4.6759
x0019 =    -3.5503
y0019 =     8.2467
z0019 =    -4.9269
 Constants:
x0001 =    -2.7880
y0001 =     9.5740
z0001 =    -0.6330
x0002 =    -2.3740
y0002 =    10.4740
z0002 =    -1.7250
x0003 =    -3.2950
y0003 =    11.6670
z0003 =    -1.6220
x0004 =    -4.5090
y0004 =    11.4480
z0004 =    -1.6140
x0005 =    -2.7060
y0005 =     9.9090
z0005 =    -3.1220
x0006 =    -2.2070
y0006 =     8.4690
z0006 =    -3.2820
x0007 =    -2.1370
y0007 =    10.8400
z0007 =    -4.2040
x0008 =    -3.0590
y0008 =     7.6140
z0008 =    -4.1730
x0020 =     2.1110
y0020 =    14.1170
z0020 =    -2.6550
x0021 =     2.5570
y0021 =    13.2420
z0021 =    -1.5050
x0022 =     1.8980
y0022 =    11.8740
z0022 =    -1.6180
x0023 =     0.6860
y0023 =    11.7740
z0023 =    -1.8300
x0024 =     2.2220
y0024 =    13.9240
z0024 =    -0.1460
x0025 =     2.5150
y0025 =    13.0990
z0025 =     1.1440
x0026 =     1.8880
y0026 =    14.0590
z0026 =     2.6100
x0027 =     2.3090
y0027 =    15.7520
z0027 =     2.1120
x0028 =     2.4559
y0028 =    15.0631
z0028 =    -2.7330
x0029 =     3.6476
y0029 =    13.1053
z0029 =    -1.5492
x0030 =     1.1488
y0030 =    14.1654
z0030 =    -0.1502
x0031 =     2.8072
y0031 =    14.8535
z0031 =    -0.0868
x0032 =     3.5984
y0032 =    12.9360
z0032 =     1.2426
x0033 =     2.0011
y0033 =    12.1280
z0033 =     1.0886
x0034 =     1.7996
y0034 =    16.4662
z0034 =     2.7756
x0035 =     1.9851
y0035 =    15.9196
z0035 =     1.0742
x0036 =     3.3972
y0036 =    15.8951
z0036 =     2.1847
x0037 =     1.5122
y0037 =    14.0632
z0037 =    -3.4666
x0038 =     2.3758
y0038 =    13.6438
z0038 =    -3.5065
x0039 =     2.7090
y0039 =    10.8300
z0039 =    -1.5060
x0040 =     2.2350
y0040 =     9.4710
z0040 =    -1.4030
x0041 =     2.9350
y0041 =     8.9880
z0041 =    -0.1290
x0042 =     4.1500
y0042 =     9.1980
z0042 =     0.0210
x0043 =     2.6820
y0043 =     8.5510
z0043 =    -2.6420
x0044 =     2.0770
y0044 =     7.1290
z0044 =    -2.5690
x0045 =     2.3610
y0045 =     9.1490
z0045 =    -3.9960
x0046 =     3.6899
y0046 =    11.0705
z0046 =    -1.4991
x0047 =     1.1365
y0047 =     9.4154
z0047 =    -1.3918
x0048 =     3.7763
y0048 =     8.4910
z0048 =    -2.5481
x0049 =     2.2783
y0049 =     6.6959
z0049 =    -1.5781
x0050 =     2.5327
y0050 =     6.4969
z0050 =    -3.3454
x0051 =     0.9905
y0051 =     7.1846
z0051 =    -2.7316
x0052 =     3.2596
y0052 =     9.6348
z0052 =    -4.4040
x0053 =     1.5587
y0053 =     9.8935
z0053 =    -3.8865
x0054 =     2.0319
y0054 =     8.3523
z0054 =    -4.6794
x0055 =    -2.8270
y0055 =     4.6290
z0055 =     4.3880
x0056 =    -2.1550
y0056 =     5.2040
z0056 =     3.2200
x0057 =    -2.8280
y0057 =     6.4460
z0057 =     2.7110
x0058 =    -4.0470
y0058 =     6.4960
z0058 =     2.6330
x0059 =    -2.1130
y0059 =     4.1600
z0059 =     2.0550
x0060 =    -3.4310
y0060 =     3.5340
z0060 =     1.8780
x0061 =    -1.5640
y0061 =     4.8060
z0061 =     0.6650
x0062 =    -3.8204
y0062 =     4.6750
z0062 =     4.5642
x0063 =    -1.1404
y0063 =     5.4720
z0063 =     3.5500
x0064 =    -1.3935
y0064 =     3.3771
z0064 =     2.3369
x0065 =    -3.9477
y0065 =     3.4880
z0065 =     2.8480
x0066 =    -3.3026
y0066 =     2.5159
z0066 =     1.4817
x0067 =    -4.0277
y0067 =     4.1303
z0067 =     1.1720
x0068 =    -0.7385
y0068 =     5.4992
z0068 =     0.8843
x0069 =    -2.3794
y0069 =     5.3517
z0069 =     0.1677
x0070 =    -1.2037
y0070 =     4.0030
z0070 =     0.0052
x0071 =    -2.0350
y0071 =     7.4390
z0071 =     2.3360
x0072 =    -2.5520
y0072 =     8.5730
z0072 =     1.5760
x0073 =    -1.9280
y0073 =     8.9650
z0073 =     0.2080
x0074 =    -0.7250
y0074 =     8.7310
z0074 =    -0.0700
x0075 =    -2.7010
y0075 =     9.8710
z0075 =     2.4620
x0076 =    -4.0470
y0076 =     9.8730
z0076 =     3.0530
x0077 =    -1.6630
y0077 =     9.9360
z0077 =     3.5500
x0078 =    -1.0711
y0078 =     7.3353
z0078 =     2.6193
x0079 =    -3.5185
y0079 =     8.1380
z0079 =     1.2815
x0080 =    -2.5517
y0080 =    10.7534
z0080 =     1.8225
x0081 =    -4.5263
y0081 =     8.8990
z0081 =     2.8754
x0082 =    -4.6491
y0082 =    10.6685
z0082 =     2.5897
x0083 =    -3.9713
y0083 =    10.0530
z0083 =     4.1355
x0084 =    -1.6043
y0084 =    10.9632
z0084 =     3.9391
x0085 =    -0.6851
y0085 =     9.6418
z0085 =     3.1411
x0086 =    -1.9421
y0086 =     9.2505
z0086 =     4.3637
x0087 =    -2.8390
y0087 =    12.9040
z0087 =    -1.4650
x0088 =    -1.6210
y0088 =    13.3180
z0088 =    -0.8240
x0089 =    -2.1700
y0089 =    13.9940
z0089 =     0.3990
x0090 =    -2.9570
y0090 =    14.8940
z0090 =     0.1490
x0091 =    -0.9050
y0091 =    14.4240
z0091 =    -1.6630
x0092 =    -1.8450
y0092 =    13.7540
z0092 =     1.5600
x0093 =    -3.4690
y0093 =    13.5845
z0093 =    -1.8652
x0094 =    -0.9079
y0094 =    12.4957
z0094 =    -0.6648
x0095 =    -0.4344
y0095 =    15.1510
z0095 =    -0.9847
x0096 =    -1.6426
y0096 =    14.9376
z0096 =    -2.2971
x0097 =    -0.1345
y0097 =    13.9612
z0097 =    -2.2971
x0098 =    -5.0140
y0098 =     8.4720
z0098 =    -9.0260
x0099 =    -4.4310
y0099 =     7.9090
z0099 =    -7.7960
x0100 =    -5.0880
y0100 =     6.6220
z0100 =    -7.3280
x0101 =    -6.2570
y0101 =     6.3860
z0101 =    -7.5750
x0102 =    -6.0046
y0102 =     8.5333
z0102 =    -9.2135
x0103 =    -4.5278
y0103 =     8.6566
z0103 =    -6.9950
x0104 =    -3.3661
y0104 =     7.7050
z0104 =    -7.9812
x0105 =    -4.3530
y0105 =     5.7780
z0105 =    -6.6340
x0106 =    -4.9800
y0106 =     4.5590
z0106 =    -6.1400
x0107 =    -4.2720
y0107 =     4.0000
z0107 =    -4.9340
x0108 =    -3.0460
y0108 =     4.0730
z0108 =    -4.8040
x0109 =    -3.3888
y0109 =     6.0435
z0109 =    -6.4924
x0110 =    -4.9649
y0110 =     3.8044
z0110 =    -6.9402
x0111 =    -6.0215
y0111 =     4.7831
z0111 =    -5.8662
x0112 =    -5.0510
y0112 =     3.4040
z0112 =    -4.0540
x0113 =    -4.4900
y0113 =     2.9010
z0113 =    -2.8160
x0114 =    -5.1280
y0114 =     1.5840
z0114 =    -2.4870
x0115 =    -6.3020
y0115 =     1.4560
z0115 =    -2.6880
x0116 =    -4.8400
y0116 =     3.8640
z0116 =    -1.7040
x0117 =    -4.2280
y0117 =     5.2080
z0117 =    -1.9950
x0118 =    -6.3730
y0118 =     4.0470
z0118 =    -1.6390
x0119 =    -6.0246
y0119 =     3.3364
z0119 =    -4.3141
x0120 =    -3.4006
y0120 =     2.7905
z0120 =    -2.9209
x0121 =    -4.4603
y0121 =     3.4623
z0121 =    -0.7530
x0122 =    -3.2200
y0122 =     5.2539
z0122 =    -1.5569
x0123 =    -4.8551
y0123 =     5.9985
z0123 =    -1.5569
x0124 =    -4.1622
y0124 =     5.3524
z0124 =    -3.0835
x0125 =    -6.7376
y0125 =     3.7332
z0125 =    -0.6497
x0126 =    -6.8485
y0126 =     3.4323
z0126 =    -2.4175
x0127 =    -6.6227
y0127 =     5.1056
z0127 =    -1.8036
x0128 =    -7.3590
y0128 =    14.1170
z0128 =    -2.6550
x0129 =    -6.9130
y0129 =    13.2420
z0129 =    -1.5050
x0130 =    -7.5720
y0130 =    11.8740
z0130 =    -1.6180
x0131 =    -8.7840
y0131 =    11.7740
z0131 =    -1.8300
x0132 =    -7.2480
y0132 =    13.9240
z0132 =    -0.1460
x0133 =    -6.9550
y0133 =    13.0990
z0133 =     1.1440
x0134 =    -7.5820
y0134 =    14.0590
z0134 =     2.6100
x0135 =    -7.1610
y0135 =    15.7520
z0135 =     2.1120
x0136 =    -7.0141
y0136 =    15.0631
z0136 =    -2.7330
x0137 =    -5.8224
y0137 =    13.1053
z0137 =    -1.5492
x0138 =    -8.3212
y0138 =    14.1654
z0138 =    -0.1502
x0139 =    -6.6628
y0139 =    14.8535
z0139 =    -0.0868
x0140 =    -5.8716
y0140 =    12.9360
z0140 =     1.2426
x0141 =    -7.4689
y0141 =    12.1280
z0141 =     1.0886
x0142 =    -7.6704
y0142 =    16.4662
z0142 =     2.7756
x0143 =    -7.4849
y0143 =    15.9196
z0143 =     1.0742
x0144 =    -6.0728
y0144 =    15.8951
z0144 =     2.1847
x0145 =    -7.9578
y0145 =    14.0632
z0145 =    -3.4666
x0146 =    -7.0942
y0146 =    13.6438
z0146 =    -3.5065
x0147 =    -6.7610
y0147 =    10.8300
z0147 =    -1.5060
x0148 =    -7.2350
y0148 =     9.4710
z0148 =    -1.4030
x0149 =    -6.5350
y0149 =     8.9880
z0149 =    -0.1290
x0150 =    -5.3200
y0150 =     9.1980
z0150 =     0.0210
x0151 =    -6.7880
y0151 =     8.5510
z0151 =    -2.6420
x0152 =    -7.3930
y0152 =     7.1290
z0152 =    -2.5690
x0153 =    -7.1090
y0153 =     9.1490
z0153 =    -3.9960
x0154 =    -5.7801
y0154 =    11.0705
z0154 =    -1.4991
x0155 =    -8.3335
y0155 =     9.4154
z0155 =    -1.3918
x0156 =    -5.6937
y0156 =     8.4910
z0156 =    -2.5481
x0157 =    -7.1917
y0157 =     6.6959
z0157 =    -1.5781
x0158 =    -6.9373
y0158 =     6.4969
z0158 =    -3.3454
x0159 =    -8.4795
y0159 =     7.1846
z0159 =    -2.7316
x0160 =    -6.2104
y0160 =     9.6348
z0160 =    -4.4040
x0161 =    -7.9113
y0161 =     9.8935
z0161 =    -3.8865
x0162 =    -7.4381
y0162 =     8.3523
z0162 =    -4.6794
x0163 =    -0.3120
y0163 =     4.4850
z0163 =    -5.5910
x0164 =     0.2720
y0164 =     4.9660
z0164 =    -6.8800
x0165 =    -0.3190
y0165 =     6.2880
z0165 =    -7.3340
x0166 =    -1.4860
y0166 =     6.5160
z0166 =    -7.1390
x0167 =     0.1030
y0167 =     3.8540
z0167 =    -8.0010
x0168 =    -1.3100
y0168 =     3.4020
z0168 =    -8.1210
x0169 =     0.6320
y0169 =     4.2970
z0169 =    -9.3920
x0170 =    -1.2856
y0170 =     4.5575
z0170 =    -5.3322
x0171 =     1.3430
y0171 =     5.1484
z0171 =    -6.7077
x0172 =     0.7244
y0172 =     3.0103
z0172 =    -7.6663
x0173 =    -1.7909
y0173 =     3.4373
z0173 =    -7.1323
x0174 =    -1.3334
y0174 =     2.3711
z0174 =    -8.5041
x0175 =    -1.8485
y0175 =     4.0640
z0175 =    -8.8151
x0176 =     1.5740
y0176 =     4.8509
z0176 =    -9.2661
x0177 =    -0.1127
y0177 =     4.9447
z0177 =    -9.8777
x0178 =     0.8090
y0178 =     3.4088
z0178 =   -10.0163
x0179 =     0.4660
y0179 =     7.1610
z0179 =    -7.9490
x0180 =    -0.1090
y0180 =     8.3710
z0180 =    -8.5850
x0181 =     0.5440
y0181 =     8.9530
z0181 =    -9.8580
x0182 =     1.7700
y0182 =     8.8700
z0182 =   -10.0760
x0183 =    -0.2690
y0183 =     9.6130
z0183 =    -7.5880
x0184 =    -1.1180
y0184 =     9.2470
z0184 =    -6.4640
x0185 =     1.0730
y0185 =    10.1430
z0185 =    -7.0880
x0186 =     1.4462
y0186 =     6.9181
z0186 =    -7.9339
x0187 =    -1.0638
y0187 =     7.9184
z0187 =    -8.8908
x0188 =    -0.7417
y0188 =    10.4263
z0188 =    -8.1581
x0189 =    -1.4659
y0189 =    10.1585
z0189 =    -5.9559
x0190 =    -0.5442
y0190 =     8.6291
z0190 =    -5.7577
x0191 =    -1.9848
y0191 =     8.6772
z0191 =    -6.8301
x0192 =     0.9647
y0192 =    11.1985
z0192 =    -6.7979
x0193 =     1.8219
y0193 =    10.0573
z0193 =    -7.8891
x0194 =     1.3977
y0194 =     9.5544
z0194 =    -6.2173
x0195 =     2.7830
y0195 =     7.8530
z0195 =     2.0440
x0196 =     1.2362
y0196 =     8.1738
z0196 =     0.4996
x0197 =     2.1980
y0197 =     8.3200
z0197 =     0.7710
x0198 =     2.2218
y0198 =     7.0393
z0198 =     2.5266
x0199 =     2.8213
y0199 =     8.7048
z0199 =     2.7389
x0200 =     3.8041
y0200 =     7.4962
z0200 =     1.8439
x0201 =    -2.8720
y0201 =     3.5230
z0201 =     6.5160
x0202 =    -0.8870
y0202 =     3.7900
z0202 =     5.2250
x0203 =    -2.3563
y0203 =     2.6635
z0203 =     6.9691
x0204 =    -2.1180
y0204 =     3.9900
z0204 =     5.3110
x0205 =    -2.9311
y0205 =     4.3410
z0205 =     7.2491
x0206 =    -3.8874
y0206 =     3.2233
z0206 =     6.2175
x0207 =    -4.9370
y0207 =     9.4920
z0207 =   -11.2240
x0208 =    -3.0060
y0208 =     8.9700
z0208 =    -9.9360
x0209 =    -4.5256
y0209 =    10.5121
z0209 =   -11.2147
x0210 =    -4.2390
y0210 =     8.9510
z0210 =    -9.9990
x0211 =    -4.6065
y0211 =     8.9314
z0211 =   -12.1109
x0212 =    -6.0354
y0212 =     9.4345
z0212 =   -11.2369
x0213 =    -4.9200
y0213 =    -0.6090
z0213 =    -1.4000
x0214 =    -3.3853
y0214 =     0.8501
z0214 =    -1.9949
x0215 =    -4.3710
y0215 =     0.6320
z0215 =    -1.9640
x0216 =    -4.5173
y0216 =    -0.7705
z0216 =    -0.3892
x0217 =    -4.6375
y0217 =    -1.4979
z0217 =    -1.9832
x0218 =    -6.0133
y0218 =    -0.4878
z0218 =    -1.4054
x0219 =    -6.6870
y0219 =     7.8530
z0219 =     2.0440
x0220 =    -8.2338
y0220 =     8.1738
z0220 =     0.4996
x0221 =    -7.2720
y0221 =     8.3200
z0221 =     0.7710
x0222 =    -7.2482
y0222 =     7.0393
z0222 =     2.5266
x0223 =    -6.6487
y0223 =     8.7048
z0223 =     2.7389
x0224 =    -5.6659
y0224 =     7.4962
z0224 =     1.8439
x0225 =    -0.2750
y0225 =     3.4470
z0225 =    -3.4040
x0226 =     1.6840
y0226 =     3.7410
z0226 =    -4.7720
x0227 =     0.1892
y0227 =     2.5517
z0227 =    -2.9646
x0228 =     0.4470
y0228 =     3.9040
z0228 =    -4.6590
x0229 =    -0.2491
y0229 =     4.2600
z0229 =    -2.6635
x0230 =    -1.3197
y0230 =     3.2162
z0230 =    -3.6595
x0231 =     0.0470
y0231 =    10.5900
z0231 =   -11.6650
x0232 =    -1.2750
y0232 =     9.2717
z0232 =   -10.5288
x0233 =    -0.3220
y0233 =     9.5750
z0233 =   -10.6700
x0234 =    -0.7971
y0234 =    11.2677
z0234 =   -11.4696
x0235 =    -0.0596
y0235 =    10.1949
z0235 =   -12.6861
x0236 =     1.0702
y0236 =    10.9883
z0236 =   -11.5980
 
HrmStr1   * * 331. 1.09
HrmBnd1   * * * 20. 90.
AmbTrs    * * * * 0 0 0 0 0.0 0.0 1.15 0.0 3.0
VDW       MG 1.09 0.25
VDW       NT 1.8240 0.17

   1-0089 0
6-31G*  
****
  90-0090 0
6-31+G* 
****
  91-0091 0
6-31G*  
****
  92-0092 0
6-31+G* 
****
  93-0236 0
6-31G*  
****
 
END
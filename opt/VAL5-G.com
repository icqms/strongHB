g09 << +++ > VAL5-G.log
%mem=1900MB
%Nprocshared=4
%Chk=/short/d63/reimers/VAL5-G
#P ONIOM(cam-b3lyp:AMBER=hardfirst)=EmbedCharge iop(4/28=5)
   opt(z-matrix,maxcycle=7,nrscale) nosym scf(conver=6) 

  lysozyme partial opt VAL5-G

  0 1  0 1  0 1
 N-N--0.416 0 x0001 y0001 z0001    
 C-CT--0.087 0 x0002 y0002 z0002    
 C-C-0.597 0 x0003 y0003 z0003    
 O-O--0.568 0 x0004 y0004 z0004    
 C-CT-0.298 0 x0005 y0005 z0005    
 C-CT--0.319 0 x0006 y0006 z0006    
 C-CT--0.319 0 x0007 y0007 z0007    
 H-H-0.272 0 x0008 y0008 z0008    
 H-H1-0.097 0 x0009 y0009 z0009    
 H-HC--0.030 0 x0010 y0010 z0010    
 H-HC-0.079 0 x0011 y0011 z0011    
 H-HC-0.079 0 x0012 y0012 z0012    
 H-HC-0.079 0 x0013 y0013 z0013    
 H-HC-0.079 0 x0014 y0014 z0014    
 H-HC-0.079 0 x0015 y0015 z0015    
 H-HC-0.079 0 x0016 y0016 z0016    
 N-N--0.416 0 x0017 y0017 z0017 M  
 C-CT--0.025 0 x0018 y0018 z0018 M  
 C-C-0.597 0 x0019 y0019 z0019 M  
 O-O--0.568 0 x0020 y0020 z0020 M  
 H-H-0.272 0 x0021 y0021 z0021 M  
 H-H1-0.070 0 x0022 y0022 z0022 M  
 H-H1-0.070 0 x0023 y0023 z0023 M  
 N-N--0.416 0 x0024 y0024 z0024 M  
 C-CT--0.025 0 x0025 y0025 z0025 M  
 C-C-0.597 0 x0026 y0026 z0026 M  
 O-O--0.568 0 x0027 y0027 z0027 M  
 H-H-0.272 0 x0028 y0028 z0028 M  
 H-H1-0.070 0 x0029 y0029 z0029 M  
 H-H1-0.070 0 x0030 y0030 z0030 M  
 N-N--0.313 0 x0031 y0031 z0031 M  
 C-CT--0.060 0 x0032 y0032 z0032 M  
 C-C-0.544 0 x0033 y0033 z0033 M  
 O-O--0.380 0 x0034 y0034 z0034 M  
 C-CT-0.130 0 x0035 y0035 z0035 M  
 C-CT--0.043 0 x0036 y0036 z0036 M  
 C-CT--0.320 0 x0037 y0037 z0037 M  
 C-CT--0.066 0 x0038 y0038 z0038 M  
 H-H-0.150 0 x0039 y0039 z0039 M  
 H-H1-0.087 0 x0040 y0040 z0040 M  
 H-HC-0.019 0 x0041 y0041 z0041 M  
 H-HC-0.024 0 x0042 y0042 z0042 M  
 H-HC-0.024 0 x0043 y0043 z0043 M  
 H-HC-0.088 0 x0044 y0044 z0044 M  
 H-HC-0.088 0 x0045 y0045 z0045 M  
 H-HC-0.088 0 x0046 y0046 z0046 M  
 H-HC-0.019 0 x0047 y0047 z0047 M  
 H-HC-0.019 0 x0048 y0048 z0048 M  
 H-HC-0.019 0 x0049 y0049 z0049 M  
 N-N--0.416 0 x0050 y0050 z0050 M H-H1--0.144 
 C-CT--0.025 0 x0051 y0051 z0051    
 C-C-0.597 0 x0052 y0052 z0052    
 O-O--0.568 0 x0053 y0053 z0053    
 H-H-0.272 0 x0054 y0054 z0054 M  
 H-H1-0.070 0 x0055 y0055 z0055    
 H-H1-0.070 0 x0056 y0056 z0056    
 N-N--0.416 0 x0057 y0057 z0057    
 C-CT--0.025 0 x0058 y0058 z0058    
 C-C-0.597 0 x0059 y0059 z0059    
 O-O--0.568 0 x0060 y0060 z0060    
 H-H-0.272 0 x0061 y0061 z0061    
 H-H1-0.070 0 x0062 y0062 z0062    
 H-H1-0.070 0 x0063 y0063 z0063    
 N-N--0.416 0 x0064 y0064 z0064    
 C-CT--0.087 0 x0065 y0065 z0065    
 C-C-0.597 0 x0066 y0066 z0066    
 O-O--0.568 0 x0067 y0067 z0067    
 C-CT-0.298 0 x0068 y0068 z0068    
 C-CT--0.319 0 x0069 y0069 z0069    
 C-CT--0.319 0 x0070 y0070 z0070    
 H-H-0.272 0 x0071 y0071 z0071    
 H-H1-0.097 0 x0072 y0072 z0072    
 H-HC--0.030 0 x0073 y0073 z0073    
 H-HC-0.079 0 x0074 y0074 z0074    
 H-HC-0.079 0 x0075 y0075 z0075    
 H-HC-0.079 0 x0076 y0076 z0076    
 H-HC-0.079 0 x0077 y0077 z0077    
 H-HC-0.079 0 x0078 y0078 z0078    
 H-HC-0.079 0 x0079 y0079 z0079    
 N-N--0.416 0 x0080 y0080 z0080    
 C-CT--0.060 0 x0081 y0081 z0081    
 C-C-0.597 0 x0082 y0082 z0082 M H-H1-0.058  
 O-O--0.568 0 x0083 y0083 z0083 M  
 C-CT-0.130 0 x0084 y0084 z0084 M H-H1-0.058  
 C-CT--0.043 0 x0085 y0085 z0085 M  
 C-CT--0.320 0 x0086 y0086 z0086 M  
 C-CT--0.066 0 x0087 y0087 z0087 M  
 H-H-0.272 0 x0088 y0088 z0088    
 H-H1-0.087 0 x0089 y0089 z0089    
 H-HC-0.019 0 x0090 y0090 z0090 M  
 H-HC-0.024 0 x0091 y0091 z0091 M  
 H-HC-0.024 0 x0092 y0092 z0092 M  
 H-HC-0.088 0 x0093 y0093 z0093 M  
 H-HC-0.088 0 x0094 y0094 z0094 M  
 H-HC-0.088 0 x0095 y0095 z0095 M  
 H-HC-0.019 0 x0096 y0096 z0096 M  
 H-HC-0.019 0 x0097 y0097 z0097 M  
 H-HC-0.019 0 x0098 y0098 z0098 M  
 N-N--0.367 0 x0099 y0099 z0099 M  
 C-CT--0.087 0 x0100 y0100 z0100 M  
 C-C-0.597 0 x0101 y0101 z0101 M  
 O-O--0.542 0 x0102 y0102 z0102 M  
 C-CT-0.298 0 x0103 y0103 z0103 M  
 C-CT--0.319 0 x0104 y0104 z0104 M  
 C-CT--0.319 0 x0105 y0105 z0105 M  
 H-H-0.151 0 x0106 y0106 z0106 M  
 H-H1-0.097 0 x0107 y0107 z0107 M  
 H-HC--0.030 0 x0108 y0108 z0108 M  
 H-HC-0.079 0 x0109 y0109 z0109 M  
 H-HC-0.079 0 x0110 y0110 z0110 M  
 H-HC-0.079 0 x0111 y0111 z0111 M  
 H-HC-0.079 0 x0112 y0112 z0112 M  
 H-HC-0.079 0 x0113 y0113 z0113 M  
 H-HC-0.079 0 x0114 y0114 z0114 M  
 N-N--0.416 0 x0115 y0115 z0115 M  
 C-CT--0.087 0 x0116 y0116 z0116 M  
 C-C-0.597 0 x0117 y0117 z0117 M  
 O-O--0.568 0 x0118 y0118 z0118 M  
 C-CT-0.169 0 x0119 y0119 z0119 M  
 C-CT--0.189 0 x0120 y0120 z0120 M  
 C-CT--0.073 0 x0121 y0121 z0121 M  
 H-H-0.272 0 x0122 y0122 z0122 M  
 H-H1-0.097 0 x0123 y0123 z0123 M  
 H-HC--0.009 0 x0124 y0124 z0124 M  
 H-HC-0.056 0 x0125 y0125 z0125 M  
 H-HC-0.074 0 x0126 y0126 z0126 M  
 H-HC-0.017 0 x0127 y0127 z0127 M  
 H-HC-0.015 0 x0128 y0128 z0128 M  
 H-HC-0.023 0 x0129 y0129 z0129 M  
 H-HC-0.009 0 x0130 y0130 z0130 M  
 N-N--0.416 0 x0131 y0131 z0131 M  
 C-CT--0.060 0 x0132 y0132 z0132 M  
 C-C-0.226 0 x0133 y0133 z0133 M  
 O-O--0.119 0 x0134 y0134 z0134 M  
 C-CT-0.130 0 x0135 y0135 z0135 M  
 C-CT--0.043 0 x0136 y0136 z0136 M  
 C-CT--0.320 0 x0137 y0137 z0137 M  
 C-CT--0.066 0 x0138 y0138 z0138 M  
 H-H-0.272 0 x0139 y0139 z0139 M  
 H-H1-0.087 0 x0140 y0140 z0140 M  
 H-HC-0.017 0 x0141 y0141 z0141 M  
 H-HC-0.024 0 x0142 y0142 z0142 M  
 H-HC-0.024 0 x0143 y0143 z0143 M  
 H-HC-0.088 0 x0144 y0144 z0144 M  
 H-HC-0.088 0 x0145 y0145 z0145 M  
 H-HC-0.088 0 x0146 y0146 z0146 M  
 H-HC-0.019 0 x0147 y0147 z0147 M  
 H-HC-0.019 0 x0148 y0148 z0148 M  
 H-HC-0.019 0 x0149 y0149 z0149 M  
 N-N--0.416 0 x0150 y0150 z0150 M  
 C-CT--0.087 0 x0151 y0151 z0151 M  
 C-C-0.597 0 x0152 y0152 z0152 M  
 O-O--0.568 0 x0153 y0153 z0153 M  
 C-CT-0.298 0 x0154 y0154 z0154 M  
 C-CT--0.319 0 x0155 y0155 z0155 M  
 C-CT--0.319 0 x0156 y0156 z0156 M  
 H-H-0.272 0 x0157 y0157 z0157 M  
 H-H1-0.097 0 x0158 y0158 z0158 M  
 H-HC--0.030 0 x0159 y0159 z0159 M  
 H-HC-0.079 0 x0160 y0160 z0160 M  
 H-HC-0.077 0 x0161 y0161 z0161 M  
 H-HC-0.077 0 x0162 y0162 z0162 M  
 H-HC-0.079 0 x0163 y0163 z0163 M  
 H-HC-0.067 0 x0164 y0164 z0164 M  
 H-HC-0.053 0 x0165 y0165 z0165 M  
 N-N--0.416 0 x0166 y0166 z0166 M H-H1--0.144 
 C-CT--0.025 0 x0167 y0167 z0167    
 C-C-0.597 0 x0168 y0168 z0168    
 O-O--0.568 0 x0169 y0169 z0169    
 H-H-0.272 0 x0170 y0170 z0170 M  
 H-H1-0.070 0 x0171 y0171 z0171    
 H-H1-0.070 0 x0172 y0172 z0172    
 N-N--0.416 0 x0173 y0173 z0173    
 C-CT--0.025 0 x0174 y0174 z0174    
 C-C-0.597 0 x0175 y0175 z0175    
 O-O--0.568 0 x0176 y0176 z0176    
 H-H-0.272 0 x0177 y0177 z0177    
 H-H1-0.070 0 x0178 y0178 z0178    
 H-H1-0.070 0 x0179 y0179 z0179    
 N-N--0.153 0 x0180 y0180 z0180 M  
 C-CT--0.087 0 x0181 y0181 z0181 M  
 C-C-0.597 0 x0182 y0182 z0182 M  
 O-O--0.435 0 x0183 y0183 z0183 M  
 C-CT-0.298 0 x0184 y0184 z0184 M  
 C-CT--0.319 0 x0185 y0185 z0185 M  
 C-CT--0.319 0 x0186 y0186 z0186 M  
 H-H-0.064 0 x0187 y0187 z0187 M  
 H-H1-0.097 0 x0188 y0188 z0188 M  
 H-HC--0.030 0 x0189 y0189 z0189 M  
 H-HC-0.079 0 x0190 y0190 z0190 M  
 H-HC-0.079 0 x0191 y0191 z0191 M  
 H-HC-0.079 0 x0192 y0192 z0192 M  
 H-HC-0.079 0 x0193 y0193 z0193 M  
 H-HC-0.079 0 x0194 y0194 z0194 M  
 H-HC-0.079 0 x0195 y0195 z0195 M  
 N-N--0.416 0 x0196 y0196 z0196 M  
 C-CT--0.025 0 x0197 y0197 z0197 M  
 C-C-0.597 0 x0198 y0198 z0198 M  
 O-O--0.568 0 x0199 y0199 z0199 M  
 H-H-0.272 0 x0200 y0200 z0200 M  
 H-H1-0.054 0 x0201 y0201 z0201 M  
 H-H1-0.070 0 x0202 y0202 z0202 M  
 N-N--0.416 0 x0203 y0203 z0203 M  
 C-CT--0.025 0 x0204 y0204 z0204 M  
 C-C-0.597 0 x0205 y0205 z0205 M  
 O-O--0.422 0 x0206 y0206 z0206 M  
 H-H-0.272 0 x0207 y0207 z0207 M  
 H-H1-0.070 0 x0208 y0208 z0208 M  
 H-H1-0.070 0 x0209 y0209 z0209 M  
 C-CT--0.087 0 x0210 y0210 z0210 M  
 O-O--0.568 0 x0211 y0211 z0211 M  
 H-N--0.416 0 x0212 y0212 z0212 M  
 C-C-0.597 0 x0213 y0213 z0213 M  
 H-CT-0.298 0 x0214 y0214 z0214 M  
 H-H1-0.097 0 x0215 y0215 z0215 M  
 C-CT--0.087 0 x0216 y0216 z0216 M  
 H-H-0.272 0 x0217 y0217 z0217 M  
 N-N--0.416 0 x0218 y0218 z0218 M  
 H-C-0.597 0 x0219 y0219 z0219 M  
 H-CT-0.298 0 x0220 y0220 z0220 M  
 H-H1-0.097 0 x0221 y0221 z0221 M  
 C-CT--0.034 0 x0222 y0222 z0222 M  
 O-O--0.497 0 x0223 y0223 z0223 M  
 H-N--0.207 0 x0224 y0224 z0224 M  
 C-C-0.387 0 x0225 y0225 z0225 M  
 H-CT-0.049 0 x0226 y0226 z0226 M  
 H-H1-0.080 0 x0227 y0227 z0227 M  
 C-CT-0.012 0 x0228 y0228 z0228 M  
 H-H-0.074 0 x0229 y0229 z0229 M  
 N-N--0.200 0 x0230 y0230 z0230 M  
 H-C-0.058 0 x0231 y0231 z0231 M  
 H-CT--0.082 0 x0232 y0232 z0232 M  
 H-H1-0.036 0 x0233 y0233 z0233 M  
 C-CT--0.087 0 x0234 y0234 z0234 M  
 O-O--0.568 0 x0235 y0235 z0235 M  
 H-N--0.416 0 x0236 y0236 z0236 M  
 C-C-0.597 0 x0237 y0237 z0237 M  
 H-CT-0.298 0 x0238 y0238 z0238 M  
 H-H1-0.097 0 x0239 y0239 z0239 M  
 C-CT-0.034 0 x0240 y0240 z0240 M  
 H-H-0.272 0 x0241 y0241 z0241 M  
 N-N--0.416 0 x0242 y0242 z0242 M  
 H-C-0.412 0 x0243 y0243 z0243 M  
 H-CT--0.182 0 x0244 y0244 z0244 M  
 H-H1-0.082 0 x0245 y0245 z0245 M  
 C-CT--0.006 0 x0246 y0246 z0246 M  
 O-O--0.555 0 x0247 y0247 z0247 M  
 H-N--0.082 0 x0248 y0248 z0248 M  
 C-C-0.437 0 x0249 y0249 z0249 M  
 H-H1-0.008 0 x0250 y0250 z0250 M  
 H-H1-0.013 0 x0251 y0251 z0251 M  
 C-CT-0.005 0 x0252 y0252 z0252 M  
 H-H-0.166 0 x0253 y0253 z0253 M  
 N-N--0.130 0 x0254 y0254 z0254 M  
 H-C-0.055 0 x0255 y0255 z0255 M  
 H-CT--0.014 0 x0256 y0256 z0256 M  
 H-H1-0.010 0 x0257 y0257 z0257 M  
 C-CT--0.023 0 x0258 y0258 z0258 M  
 O-O--0.363 0 x0259 y0259 z0259 M  
 H-N3-0.090 0 x0260 y0260 z0260 M  
 C-C-0.597 0 x0261 y0261 z0261 M  
 H-CT-0.034 0 x0262 y0262 z0262 M  
 H-H1-0.088 0 x0263 y0263 z0263 M  
 C-CT--0.087 0 x0264 y0264 z0264    
 H-H-0.272 0 x0265 y0265 z0265    
 N-N--0.416 0 x0266 y0266 z0266    
 H-C-0.597 0 x0267 y0267 z0267    
 H-CT-0.298 0 x0268 y0268 z0268    
 H-H1-0.097 0 x0269 y0269 z0269    
 C-CT--0.005 0 x0270 y0270 z0270 M  
 O-O--0.212 0 x0271 y0271 z0271 M  
 H-N3-0.061 0 x0272 y0272 z0272 M  
 C-C-0.205 0 x0273 y0273 z0273 M  
 H-CT-0.004 0 x0274 y0274 z0274 M  
 H-H1-0.015 0 x0275 y0275 z0275 M  
 C-CT--0.016 0 x0276 y0276 z0276 M  
 H-H-0.240 0 x0277 y0277 z0277 M  
 N-N--0.282 0 x0278 y0278 z0278 M  
 H-C-0.088 0 x0279 y0279 z0279 M  
 H-CT-0.034 0 x0280 y0280 z0280 M  
 H-H1-0.016 0 x0281 y0281 z0281 M  
 Variables:
x0008 =    -3.8204
y0008 =     4.6750
z0008 =     4.5642
x0009 =    -1.1404
y0009 =     5.4720
z0009 =     3.5500
x0010 =    -1.3935
y0010 =     3.3771
z0010 =     2.3369
x0011 =    -3.9477
y0011 =     3.4880
z0011 =     2.8480
x0012 =    -3.3026
y0012 =     2.5159
z0012 =     1.4817
x0013 =    -4.0277
y0013 =     4.1303
z0013 =     1.1720
x0014 =    -0.7385
y0014 =     5.4992
z0014 =     0.8843
x0015 =    -2.3794
y0015 =     5.3517
z0015 =     0.1677
x0016 =    -1.2037
y0016 =     4.0030
z0016 =     0.0052
 Constants:
x0001 =    -2.8270
y0001 =     4.6290
z0001 =     4.3880
x0002 =    -2.1550
y0002 =     5.2040
z0002 =     3.2200
x0003 =    -2.8280
y0003 =     6.4460
z0003 =     2.7110
x0004 =    -4.0470
y0004 =     6.4960
z0004 =     2.6330
x0005 =    -2.1130
y0005 =     4.1600
z0005 =     2.0550
x0006 =    -3.4310
y0006 =     3.5340
z0006 =     1.8780
x0007 =    -1.5640
y0007 =     4.8060
z0007 =     0.6650
x0017 =     2.1980
y0017 =     8.3200
z0017 =     0.7710
x0018 =     2.7830
y0018 =     7.8530
z0018 =     2.0440
x0019 =     2.0050
y0019 =     6.7250
z0019 =     2.7130
x0020 =     0.7990
y0020 =     6.6320
z0020 =     2.5430
x0021 =     1.2362
y0021 =     8.1738
z0021 =     0.4996
x0022 =     2.8213
y0022 =     8.7048
z0022 =     2.7389
x0023 =     3.8041
y0023 =     7.4962
z0023 =     1.8439
x0024 =     2.6790
y0024 =     5.8770
z0024 =     3.4880
x0025 =     2.0050
y0025 =     4.6830
z0025 =     4.0380
x0026 =     2.7060
y0026 =     4.0060
z0026 =     5.2120
x0027 =     3.9310
y0027 =     4.0070
z0027 =     5.3390
x0028 =     3.6447
y0028 =     6.1249
z0028 =     3.6492
x0029 =     1.9136
y0029 =     3.9442
z0029 =     3.2282
x0030 =     1.0022
y0030 =     4.9862
z0030 =     4.3733
x0031 =    -5.2020
y0031 =     9.5750
z0031 =     9.5490
x0032 =    -4.8330
y0032 =    10.5900
z0032 =     8.5540
x0033 =    -5.9950
y0033 =    11.5230
z0033 =     8.8230
x0034 =    -7.1170
y0034 =    11.0270
z0034 =     8.7910
x0035 =    -4.9830
y0035 =    10.0340
z0035 =     7.1170
x0036 =    -4.5840
y0036 =     8.5530
z0036 =     7.0580
x0037 =    -4.2140
y0037 =    10.9200
z0037 =     6.0910
x0038 =    -5.0540
y0038 =     7.7620
z0038 =     5.8560
x0039 =    -6.1550
y0039 =     9.2717
z0039 =     9.6902
x0040 =    -3.8098
y0040 =    10.9883
z0040 =     8.6210
x0041 =    -6.0434
y0041 =    10.0804
z0041 =     6.8283
x0042 =    -4.9896
y0042 =     8.0641
z0042 =     7.9560
x0043 =    -3.4853
y0043 =     8.5044
z0043 =     7.0781
x0044 =    -3.4565
y0044 =    11.5179
z0044 =     6.6189
x0045 =    -4.9225
y0045 =    11.5907
z0045 =     5.5829
x0046 =    -3.7210
y0046 =    10.2758
z0046 =     5.3481
x0047 =    -4.4736
y0047 =     6.8306
z0047 =     5.7809
x0048 =    -4.9082
y0048 =     8.3600
z0048 =     4.9443
x0049 =    -6.1212
y0049 =     7.5214
z0049 =     5.9705
x0050 =    -2.7020
y0050 =    -0.3380
z0050 =     9.4720
x0051 =    -2.1860
y0051 =     0.2550
z0051 =     8.2550
x0052 =    -2.8650
y0052 =     1.5690
z0052 =     7.9020
x0053 =    -3.9910
y0053 =     1.8650
z0053 =     8.3320
x0054 =    -3.6923
y0054 =    -0.4051
z0054 =     9.6587
x0055 =    -2.3401
y0055 =    -0.4529
z0055 =     7.4272
x0056 =    -1.1094
y0056 =     0.4402
z0056 =     8.3842
x0057 =    -2.1880
y0057 =     2.3830
z0057 =     7.1170
x0058 =    -2.8720
y0058 =     3.5230
z0058 =     6.5160
x0059 =    -2.1180
y0059 =     3.9900
z0059 =     5.3110
x0060 =    -0.8870
y0060 =     3.7900
z0060 =     5.2250
x0061 =    -1.2142
y0061 =     2.1497
z0061 =     6.9853
x0062 =    -2.9311
y0062 =     4.3410
z0062 =     7.2491
x0063 =    -3.8874
y0063 =     3.2233
z0063 =     6.2175
x0064 =    -2.0350
y0064 =     7.4390
z0064 =     2.3360
x0065 =    -2.5520
y0065 =     8.5730
z0065 =     1.5760
x0066 =    -1.9280
y0066 =     8.9650
z0066 =     0.2080
x0067 =    -0.7250
y0067 =     8.7310
z0067 =    -0.0700
x0068 =    -2.7010
y0068 =     9.8710
z0068 =     2.4620
x0069 =    -4.0470
y0069 =     9.8730
z0069 =     3.0530
x0070 =    -1.6630
y0070 =     9.9360
z0070 =     3.5500
x0071 =    -1.0711
y0071 =     7.3353
z0071 =     2.6193
x0072 =    -3.5185
y0072 =     8.1380
z0072 =     1.2815
x0073 =    -2.5517
y0073 =    10.7534
z0073 =     1.8225
x0074 =    -4.5263
y0074 =     8.8990
z0074 =     2.8754
x0075 =    -4.6491
y0075 =    10.6685
z0075 =     2.5897
x0076 =    -3.9713
y0076 =    10.0530
z0076 =     4.1355
x0077 =    -1.6043
y0077 =    10.9632
z0077 =     3.9391
x0078 =    -0.6851
y0078 =     9.6418
z0078 =     3.1411
x0079 =    -1.9421
y0079 =     9.2505
z0079 =     4.3637
x0080 =    -2.7880
y0080 =     9.5740
z0080 =    -0.6330
x0081 =    -2.3740
y0081 =    10.4740
z0081 =    -1.7250
x0082 =    -3.2950
y0082 =    11.6670
z0082 =    -1.6220
x0083 =    -4.5090
y0083 =    11.4480
z0083 =    -1.6140
x0084 =    -2.7060
y0084 =     9.9090
z0084 =    -3.1220
x0085 =    -2.2070
y0085 =     8.4690
z0085 =    -3.2820
x0086 =    -2.1370
y0086 =    10.8400
z0086 =    -4.2040
x0087 =    -3.0590
y0087 =     7.6140
z0087 =    -4.1730
x0088 =    -3.7550
y0088 =     9.3531
z0088 =    -0.4425
x0089 =    -1.2934
y0089 =    10.6566
z0089 =    -1.6307
x0090 =    -3.7992
y0090 =     9.8702
z0090 =    -3.2381
x0091 =    -2.1753
y0091 =     8.0032
z0091 =    -2.2860
x0092 =    -1.1919
y0092 =     8.5016
z0092 =    -3.7046
x0093 =    -2.7028
y0093 =    11.7833
z0093 =    -4.2099
x0094 =    -2.2228
y0094 =    10.3542
z0094 =    -5.1872
x0095 =    -1.0786
y0095 =    11.0481
z0095 =    -3.9884
x0096 =    -3.8227
y0096 =     7.1027
z0096 =    -3.5686
x0097 =    -2.4278
y0097 =     6.8666
z0097 =    -4.6759
x0098 =    -3.5503
y0098 =     8.2467
z0098 =    -4.9269
x0099 =    -5.0510
y0099 =     3.4040
z0099 =    -4.0540
x0100 =    -4.4900
y0100 =     2.9010
z0100 =    -2.8160
x0101 =    -5.1280
y0101 =     1.5840
z0101 =    -2.4870
x0102 =    -6.3020
y0102 =     1.4560
z0102 =    -2.6880
x0103 =    -4.8400
y0103 =     3.8640
z0103 =    -1.7040
x0104 =    -4.2280
y0104 =     5.2080
z0104 =    -1.9950
x0105 =    -6.3730
y0105 =     4.0470
z0105 =    -1.6390
x0106 =    -6.0246
y0106 =     3.3364
z0106 =    -4.3141
x0107 =    -3.4006
y0107 =     2.7905
z0107 =    -2.9209
x0108 =    -4.4603
y0108 =     3.4623
z0108 =    -0.7530
x0109 =    -3.2200
y0109 =     5.2539
z0109 =    -1.5569
x0110 =    -4.8551
y0110 =     5.9985
z0110 =    -1.5569
x0111 =    -4.1622
y0111 =     5.3524
z0111 =    -3.0835
x0112 =    -6.7376
y0112 =     3.7332
z0112 =    -0.6497
x0113 =    -6.8485
y0113 =     3.4323
z0113 =    -2.4175
x0114 =    -6.6227
y0114 =     5.1056
z0114 =    -1.8036
x0115 =    -4.3710
y0115 =     0.6320
z0115 =    -1.9640
x0116 =    -4.9200
y0116 =    -0.6090
z0116 =    -1.4000
x0117 =    -4.3690
y0117 =    -0.8300
z0117 =    -0.0170
x0118 =    -3.2060
y0118 =    -0.4990
z0118 =     0.2790
x0119 =    -4.5180
y0119 =    -1.8740
z0119 =    -2.2300
x0120 =    -3.1140
y0120 =    -2.1950
z0120 =    -1.9450
x0121 =    -5.3030
y0121 =    -3.0980
z0121 =    -1.8520
x0122 =    -3.3853
y0122 =     0.8501
z0122 =    -1.9949
x0123 =    -6.0133
y0123 =    -0.4878
z0123 =    -1.4054
x0124 =    -4.7095
y0124 =    -1.6304
z0124 =    -3.2855
x0125 =    -3.0274
y0125 =    -2.5787
z0125 =    -0.9177
x0126 =    -2.5021
y0126 =    -1.2870
z0126 =    -2.0506
x0127 =    -2.7621
y0127 =    -2.9594
z0127 =    -2.6534
x0128 =    -6.3786
y0128 =    -2.8830
z0128 =    -1.9342
x0129 =    -5.0633
y0129 =    -3.3814
z0129 =    -0.8165
x0130 =    -5.0410
y0130 =    -3.9245
z0130 =    -2.5289
x0131 =    -5.2070
y0131 =    -1.4240
z0131 =     0.8350
x0132 =    -4.7940
y0132 =    -1.9670
z0132 =     2.1570
x0133 =    -5.3050
y0133 =    -3.3900
z0133 =     2.2420
x0134 =    -6.4490
y0134 =    -3.6610
z0134 =     1.9120
x0135 =    -5.4110
y0135 =    -1.2060
z0135 =     3.3400
x0136 =    -5.0910
y0136 =     0.2810
z0136 =     3.2860
x0137 =    -4.9020
y0137 =    -1.8050
z0137 =     4.6740
x0138 =    -5.5310
y0138 =     1.0210
z0138 =     4.5040
x0139 =    -6.1518
y0139 =    -1.4625
z0139 =     0.4800
x0140 =    -3.6994
y0140 =    -1.8817
z0140 =     2.2241
x0141 =    -6.5036
y0141 =    -1.3160
z0141 =     3.2751
x0142 =    -5.5966
y0142 =     0.7158
z0142 =     2.4112
x0143 =    -4.0028
y0143 =     0.4007
z0143 =     3.1783
x0144 =    -4.1037
y0144 =    -1.1669
z0144 =     5.0809
x0145 =    -4.5082
y0145 =    -2.8163
z0145 =     4.4946
x0146 =    -5.7330
y0146 =    -1.8567
z0146 =     5.3928
x0147 =    -5.0579
y0147 =     0.5777
z0147 =     5.3926
x0148 =    -6.6249
y0148 =     0.9535
z0148 =     4.5982
x0149 =    -5.2343
y0149 =     2.0767
z0149 =     4.4181
x0150 =    -6.7610
y0150 =    10.8300
z0150 =    -1.5060
x0151 =    -7.2350
y0151 =     9.4710
z0151 =    -1.4030
x0152 =    -6.5350
y0152 =     8.9880
z0152 =    -0.1290
x0153 =    -5.3200
y0153 =     9.1980
z0153 =     0.0210
x0154 =    -6.7880
y0154 =     8.5510
z0154 =    -2.6420
x0155 =    -7.3930
y0155 =     7.1290
z0155 =    -2.5690
x0156 =    -7.1090
y0156 =     9.1490
z0156 =    -3.9960
x0157 =    -5.7801
y0157 =    11.0705
z0157 =    -1.4991
x0158 =    -8.3335
y0158 =     9.4154
z0158 =    -1.3918
x0159 =    -5.6937
y0159 =     8.4910
z0159 =    -2.5481
x0160 =    -7.1917
y0160 =     6.6959
z0160 =    -1.5781
x0161 =    -6.9373
y0161 =     6.4969
z0161 =    -3.3454
x0162 =    -8.4795
y0162 =     7.1846
z0162 =    -2.7316
x0163 =    -6.2104
y0163 =     9.6348
z0163 =    -4.4040
x0164 =    -7.9113
y0164 =     9.8935
z0164 =    -3.8865
x0165 =    -7.4381
y0165 =     8.3523
z0165 =    -4.6794
x0166 =    -7.2720
y0166 =     8.3200
z0166 =     0.7710
x0167 =    -6.6870
y0167 =     7.8530
z0167 =     2.0440
x0168 =    -7.4650
y0168 =     6.7250
z0168 =     2.7130
x0169 =    -8.6710
y0169 =     6.6320
z0169 =     2.5430
x0170 =    -8.2338
y0170 =     8.1738
z0170 =     0.4996
x0171 =    -6.6487
y0171 =     8.7048
z0171 =     2.7389
x0172 =    -5.6659
y0172 =     7.4962
z0172 =     1.8439
x0173 =    -6.7910
y0173 =     5.8770
z0173 =     3.4880
x0174 =    -7.4650
y0174 =     4.6830
z0174 =     4.0380
x0175 =    -6.7640
y0175 =     4.0060
z0175 =     5.2120
x0176 =    -5.5390
y0176 =     4.0070
z0176 =     5.3390
x0177 =    -5.8253
y0177 =     6.1249
z0177 =     3.6492
x0178 =    -7.5564
y0178 =     3.9442
z0178 =     3.2282
x0179 =    -8.4678
y0179 =     4.9862
z0179 =     4.3733
x0180 =     0.3030
y0180 =    -2.6960
z0180 =     2.0750
x0181 =    -0.2530
y0181 =    -1.4370
z0181 =     1.7140
x0182 =     0.4990
y0182 =    -0.8860
z0182 =     0.4540
x0183 =     1.7230
y0183 =    -0.8090
z0183 =     0.4420
x0184 =    -0.0860
y0184 =    -0.4400
z0184 =     2.9440
x0185 =    -0.5620
y0185 =     1.0120
z0185 =     2.6400
x0186 =    -0.6890
y0186 =    -0.9980
z0186 =     4.1960
x0187 =     1.2898
y0187 =    -2.8727
z0187 =     2.1983
x0188 =    -1.3209
y0188 =    -1.5387
z0188 =     1.4705
x0189 =     0.9963
y0189 =    -0.3510
z0189 =     3.1192
x0190 =    -1.3557
y0190 =     1.2929
z0190 =     3.3479
x0191 =    -0.9517
y0191 =     1.0620
z0191 =     1.6126
x0192 =     0.2859
y0192 =     1.7049
z0192 =     2.7452
x0193 =     0.0433
y0193 =    -1.6484
z0193 =     4.6968
x0194 =    -1.5856
y0194 =    -1.5823
z0194 =     3.9417
x0195 =    -0.9670
y0195 =    -0.1726
z0195 =     4.8679
x0196 =    -0.2500
y0196 =    -0.5090
z0196 =    -0.5840
x0197 =     0.3030
y0197 =     0.0730
z0197 =    -1.7930
x0198 =    -0.4050
y0198 =     1.3320
z0198 =    -2.2320
x0199 =    -1.6170
y0199 =     1.4790
z0199 =    -2.0580
x0200 =    -1.2394
y0200 =    -0.6648
z0200 =    -0.4541
x0201 =     0.2287
y0201 =    -0.6685
z0201 =    -2.6021
x0202 =     1.3609
y0202 =     0.3136
z0202 =    -1.6114
x0203 =     0.3420
y0203 =     2.2570
z0203 =    -2.8200
x0204 =    -0.2750
y0204 =     3.4470
z0204 =    -3.4040
x0205 =     0.4470
y0205 =     3.9040
z0205 =    -4.6590
x0206 =     1.6840
y0206 =     3.7410
z0206 =    -4.7720
x0207 =     1.3335
y0207 =     2.0644
z0207 =    -2.8229
x0208 =    -0.2491
y0208 =     4.2600
z0208 =    -2.6635
x0209 =    -1.3197
y0209 =     3.2162
z0209 =    -3.6595
x0210 =     2.2350
y0210 =     9.4710
z0210 =    -1.4030
x0211 =     4.1500
y0211 =     9.1980
z0211 =     0.0210
x0212 =     2.5963
y0212 =    10.5070
z0212 =    -1.4815
x0213 =     2.9350
y0213 =     8.9880
z0213 =    -0.1290
x0214 =     2.5410
y0214 =     8.8411
z0214 =    -2.2513
x0215 =     1.1365
y0215 =     9.4154
z0215 =    -1.3918
x0216 =     2.4510
y0216 =     2.9370
z0216 =     7.3480
x0217 =     0.9499
y0217 =     3.3181
z0217 =     5.8014
x0218 =     1.9150
y0218 =     3.4090
z0218 =     6.0850
x0219 =     2.0013
y0219 =     1.9653
z0219 =     7.6002
x0220 =     2.2488
y0220 =     3.6770
z0220 =     8.1364
x0221 =     3.5400
y0221 =     2.8072
z0221 =     7.2624
x0222 =    -4.9890
y0222 =     8.3710
z0222 =    11.6340
x0223 =    -3.1100
y0223 =     8.8700
z0223 =    10.1430
x0224 =    -4.5625
y0224 =     7.4735
z0224 =    12.1058
x0225 =    -4.3360
y0225 =     8.9530
z0225 =    10.3610
x0226 =    -5.0990
y0226 =     9.2245
z0226 =    12.3191
x0227 =    -5.9438
y0227 =     7.9184
z0227 =    11.3282
x0228 =    -4.5880
y0228 =    13.4940
z0228 =     9.4650
x0229 =    -6.7204
y0229 =    13.2918
z0229 =     9.2306
x0230 =    -5.8360
y0230 =    12.8070
z0230 =     9.1770
x0231 =    -4.7656
y0231 =    14.1354
z0231 =    10.3408
x0232 =    -4.3298
y0232 =    14.1639
z0232 =     8.6316
x0233 =    -3.7858
y0233 =    12.7592
z0233 =     9.6280
x0234 =    -2.6900
y0234 =    -1.2810
z0234 =    11.7490
x0235 =    -0.7010
y0235 =    -1.0090
z0235 =    10.3210
x0236 =    -2.2774
y0236 =    -2.2584
z0236 =    12.0395
x0237 =    -1.9340
y0237 =    -0.8440
z0237 =    10.4470
x0238 =    -2.5321
y0238 =    -0.5744
z0238 =    12.5771
x0239 =    -3.7709
y0239 =    -1.3162
z0239 =    11.5477
x0240 =    -1.6210
y0240 =    13.3180
z0240 =    -0.8240
x0241 =    -3.4690
y0241 =    13.5845
z0241 =    -1.8652
x0242 =    -2.8390
y0242 =    12.9040
z0242 =    -1.4650
x0243 =    -2.0232
y0243 =    13.8133
z0243 =     0.0720
x0244 =    -1.1168
y0244 =    14.0969
z0244 =    -1.4148
x0245 =    -0.9079
y0245 =    12.4957
z0245 =    -0.6648
x0246 =    -4.9800
y0246 =     4.5590
z0246 =    -6.1400
x0247 =    -3.0460
y0247 =     4.0730
z0247 =    -4.8040
x0248 =    -4.5067
y0248 =     5.4793
z0248 =    -6.5129
x0249 =    -4.2720
y0249 =     4.0000
z0249 =    -4.9340
x0250 =    -4.9649
y0250 =     3.8044
z0250 =    -6.9402
x0251 =    -6.0215
y0251 =     4.7831
z0251 =    -5.8662
x0252 =    -4.9640
y0252 =    -5.6570
z0252 =     2.6390
x0253 =    -3.5947
y0253 =    -4.0038
z0253 =     3.0205
x0254 =    -4.4960
y0254 =    -4.3080
z0254 =     2.6810
x0255 =    -4.2980
y0255 =    -6.3564
z0255 =     3.1656
x0256 =    -5.0685
y0256 =    -5.9849
z0256 =     1.5942
x0257 =    -5.9371
y0257 =    -5.6691
z0257 =     3.1518
x0258 =    -6.9130
y0258 =    13.2420
z0258 =    -1.5050
x0259 =    -8.7840
y0259 =    11.7740
z0259 =    -1.8300
x0260 =    -7.2374
y0260 =    13.8784
z0260 =    -2.3415
x0261 =    -7.5720
y0261 =    11.8740
z0261 =    -1.6180
x0262 =    -7.1497
y0262 =    13.7238
z0262 =    -0.5449
x0263 =    -5.8224
y0263 =    13.1053
z0263 =    -1.5492
x0264 =    -7.0190
y0264 =     2.9370
z0264 =     7.3480
x0265 =    -8.5201
y0265 =     3.3181
z0265 =     5.8014
x0266 =    -7.5550
y0266 =     3.4090
z0266 =     6.0850
x0267 =    -7.4687
y0267 =     1.9653
z0267 =     7.6002
x0268 =    -7.2212
y0268 =     3.6770
z0268 =     8.1364
x0269 =    -5.9300
y0269 =     2.8072
z0269 =     7.2624
x0270 =     0.2040
y0270 =    -4.9640
z0270 =     2.8010
x0271 =    -1.6900
y0271 =    -3.6770
z0271 =     2.0570
x0272 =    -0.3459
y0272 =    -5.2255
z0272 =     3.7171
x0273 =    -0.4700
y0273 =    -3.7170
z0273 =     2.2700
x0274 =     0.1277
y0274 =    -5.8010
z0274 =     2.0914
x0275 =     1.2747
y0275 =    -4.7795
z0275 =     2.9732
x0276 =     0.2720
y0276 =     4.9660
z0276 =    -6.8800
x0277 =    -1.2856
y0277 =     4.5575
z0277 =    -5.3322
x0278 =    -0.3120
y0278 =     4.4850
z0278 =    -5.5910
x0279 =    -0.1564
y0279 =     5.9242
z0279 =    -7.2091
x0280 =     0.1549
y0280 =     4.1957
z0280 =    -7.6565
x0281 =     1.3430
y0281 =     5.1484
z0281 =    -6.7077
 
HrmStr1   * * 331. 1.09
HrmBnd1   * * * 20. 90.
AmbTrs    * * * * 0 0 0 0 0.0 0.0 1.15 0.0 3.0
VDW       MG 1.09 0.25
VDW       NT 1.8240 0.17

   1-0281 0
6-31G*  
****
 
END
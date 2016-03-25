g09 << +++ > VAL5-A.log
%mem=1900MB
%Nprocshared=4
%Chk=/short/d63/reimers/VAL5-A
#P ONIOM(cam-b3lyp:AMBER=hardfirst)=EmbedCharge iop(4/28=5)
   opt(z-matrix,maxcycle=7,nrscale) nosym scf(conver=6) 

  lysozyme partial opt VAL5-A

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
 N-N--0.416 0 x0017 y0017 z0017 M H-H1--0.144 
 C-CT--0.025 0 x0018 y0018 z0018    
 C-C-0.597 0 x0019 y0019 z0019    
 O-O--0.568 0 x0020 y0020 z0020    
 H-H-0.272 0 x0021 y0021 z0021 M  
 H-H1-0.070 0 x0022 y0022 z0022    
 H-H1-0.070 0 x0023 y0023 z0023    
 N-N--0.416 0 x0024 y0024 z0024    
 C-CT--0.025 0 x0025 y0025 z0025    
 C-C-0.597 0 x0026 y0026 z0026    
 O-O--0.568 0 x0027 y0027 z0027    
 H-H-0.272 0 x0028 y0028 z0028    
 H-H1-0.070 0 x0029 y0029 z0029    
 H-H1-0.070 0 x0030 y0030 z0030    
 N-N--0.416 0 x0031 y0031 z0031    
 C-CT--0.087 0 x0032 y0032 z0032    
 C-C-0.597 0 x0033 y0033 z0033    
 O-O--0.568 0 x0034 y0034 z0034    
 C-CT-0.298 0 x0035 y0035 z0035    
 C-CT--0.319 0 x0036 y0036 z0036    
 C-CT--0.319 0 x0037 y0037 z0037    
 H-H-0.272 0 x0038 y0038 z0038    
 H-H1-0.097 0 x0039 y0039 z0039    
 H-HC--0.030 0 x0040 y0040 z0040    
 H-HC-0.079 0 x0041 y0041 z0041    
 H-HC-0.079 0 x0042 y0042 z0042    
 H-HC-0.079 0 x0043 y0043 z0043    
 H-HC-0.079 0 x0044 y0044 z0044    
 H-HC-0.079 0 x0045 y0045 z0045    
 H-HC-0.079 0 x0046 y0046 z0046    
 N-N--0.416 0 x0047 y0047 z0047    
 C-CT--0.060 0 x0048 y0048 z0048    
 C-C-0.597 0 x0049 y0049 z0049 M H-H1-0.058  
 O-O--0.568 0 x0050 y0050 z0050 M  
 C-CT-0.130 0 x0051 y0051 z0051 M H-H1-0.058  
 C-CT--0.043 0 x0052 y0052 z0052 M  
 C-CT--0.320 0 x0053 y0053 z0053 M  
 C-CT--0.066 0 x0054 y0054 z0054 M  
 H-H-0.272 0 x0055 y0055 z0055    
 H-H1-0.087 0 x0056 y0056 z0056    
 H-HC-0.019 0 x0057 y0057 z0057 M  
 H-HC-0.024 0 x0058 y0058 z0058 M  
 H-HC-0.024 0 x0059 y0059 z0059 M  
 H-HC-0.088 0 x0060 y0060 z0060 M  
 H-HC-0.088 0 x0061 y0061 z0061 M  
 H-HC-0.088 0 x0062 y0062 z0062 M  
 H-HC-0.019 0 x0063 y0063 z0063 M  
 H-HC-0.019 0 x0064 y0064 z0064 M  
 H-HC-0.019 0 x0065 y0065 z0065 M  
 N-N--0.195 0 x0066 y0066 z0066 M  
 C-CT--0.058 0 x0067 y0067 z0067 M  
 C-C-0.587 0 x0068 y0068 z0068 M  
 O-O--0.568 0 x0069 y0069 z0069 M  
 C-CT-0.130 0 x0070 y0070 z0070 M  
 C-CT--0.043 0 x0071 y0071 z0071 M  
 C-CT--0.320 0 x0072 y0072 z0072 M  
 C-CT--0.066 0 x0073 y0073 z0073 M  
 H-H-0.164 0 x0074 y0074 z0074 M  
 H-H1-0.072 0 x0075 y0075 z0075 M  
 H-HC-0.019 0 x0076 y0076 z0076 M  
 H-HC-0.024 0 x0077 y0077 z0077 M  
 H-HC-0.024 0 x0078 y0078 z0078 M  
 H-HC-0.088 0 x0079 y0079 z0079 M  
 H-HC-0.088 0 x0080 y0080 z0080 M  
 H-HC-0.088 0 x0081 y0081 z0081 M  
 H-HC-0.019 0 x0082 y0082 z0082 M  
 H-HC-0.019 0 x0083 y0083 z0083 M  
 H-HC-0.019 0 x0084 y0084 z0084 M  
 N-N--0.416 0 x0085 y0085 z0085 M  
 C-CT--0.087 0 x0086 y0086 z0086 M  
 C-C-0.597 0 x0087 y0087 z0087 M  
 O-O--0.568 0 x0088 y0088 z0088 M  
 C-CT-0.298 0 x0089 y0089 z0089 M  
 C-CT--0.319 0 x0090 y0090 z0090 M  
 C-CT--0.258 0 x0091 y0091 z0091 M  
 H-H-0.272 0 x0092 y0092 z0092 M  
 H-H1-0.097 0 x0093 y0093 z0093 M  
 H-HC--0.030 0 x0094 y0094 z0094 M  
 H-HC-0.067 0 x0095 y0095 z0095 M  
 H-HC-0.079 0 x0096 y0096 z0096 M  
 H-HC-0.079 0 x0097 y0097 z0097 M  
 H-HC-0.073 0 x0098 y0098 z0098 M  
 H-HC-0.039 0 x0099 y0099 z0099 M  
 H-HC-0.030 0 x0100 y0100 z0100 M  
 N-N--0.416 0 x0101 y0101 z0101 M H-H1--0.144 
 C-CT--0.025 0 x0102 y0102 z0102    
 C-C-0.597 0 x0103 y0103 z0103    
 O-O--0.568 0 x0104 y0104 z0104    
 H-H-0.272 0 x0105 y0105 z0105 M  
 H-H1-0.070 0 x0106 y0106 z0106    
 H-H1-0.070 0 x0107 y0107 z0107    
 N-N--0.416 0 x0108 y0108 z0108    
 C-CT--0.025 0 x0109 y0109 z0109    
 C-C-0.597 0 x0110 y0110 z0110    
 O-O--0.568 0 x0111 y0111 z0111    
 H-H-0.272 0 x0112 y0112 z0112    
 H-H1-0.070 0 x0113 y0113 z0113    
 H-H1-0.070 0 x0114 y0114 z0114    
 N-N--0.381 0 x0115 y0115 z0115 M  
 C-CT--0.054 0 x0116 y0116 z0116 M  
 C-C-0.415 0 x0117 y0117 z0117 M  
 O-O--0.396 0 x0118 y0118 z0118 M  
 C-CT-0.298 0 x0119 y0119 z0119 M  
 C-CT--0.319 0 x0120 y0120 z0120 M  
 C-CT--0.319 0 x0121 y0121 z0121 M  
 H-H-0.272 0 x0122 y0122 z0122 M  
 H-H1-0.066 0 x0123 y0123 z0123 M  
 H-HC--0.030 0 x0124 y0124 z0124 M  
 H-HC-0.079 0 x0125 y0125 z0125 M  
 H-HC-0.079 0 x0126 y0126 z0126 M  
 H-HC-0.079 0 x0127 y0127 z0127 M  
 H-HC-0.079 0 x0128 y0128 z0128 M  
 H-HC-0.079 0 x0129 y0129 z0129 M  
 H-HC-0.079 0 x0130 y0130 z0130 M  
 N-N--0.401 0 x0131 y0131 z0131 M  
 C-CT--0.087 0 x0132 y0132 z0132 M  
 C-C-0.597 0 x0133 y0133 z0133 M  
 O-O--0.568 0 x0134 y0134 z0134 M  
 C-CT-0.298 0 x0135 y0135 z0135 M  
 C-CT--0.319 0 x0136 y0136 z0136 M  
 C-CT--0.319 0 x0137 y0137 z0137 M  
 H-H-0.265 0 x0138 y0138 z0138 M  
 H-H1-0.097 0 x0139 y0139 z0139 M  
 H-HC--0.030 0 x0140 y0140 z0140 M  
 H-HC-0.079 0 x0141 y0141 z0141 M  
 H-HC-0.079 0 x0142 y0142 z0142 M  
 H-HC-0.079 0 x0143 y0143 z0143 M  
 H-HC-0.079 0 x0144 y0144 z0144 M  
 H-HC-0.079 0 x0145 y0145 z0145 M  
 H-HC-0.079 0 x0146 y0146 z0146 M  
 N-N--0.416 0 x0147 y0147 z0147 M  
 C-CT--0.025 0 x0148 y0148 z0148 M  
 C-C-0.597 0 x0149 y0149 z0149 M  
 O-O--0.568 0 x0150 y0150 z0150 M  
 H-H-0.272 0 x0151 y0151 z0151 M  
 H-H1-0.070 0 x0152 y0152 z0152 M  
 H-H1-0.070 0 x0153 y0153 z0153 M  
 N-N--0.416 0 x0154 y0154 z0154 M  
 C-CT--0.025 0 x0155 y0155 z0155 M  
 C-C-0.597 0 x0156 y0156 z0156 M  
 O-O--0.568 0 x0157 y0157 z0157 M  
 H-H-0.272 0 x0158 y0158 z0158 M  
 H-H1-0.070 0 x0159 y0159 z0159 M  
 H-H1-0.070 0 x0160 y0160 z0160 M  
 N-N--0.328 0 x0161 y0161 z0161 M  
 C-CT--0.087 0 x0162 y0162 z0162 M  
 C-C-0.597 0 x0163 y0163 z0163 M  
 O-O--0.568 0 x0164 y0164 z0164 M  
 C-CT-0.298 0 x0165 y0165 z0165 M  
 C-CT--0.319 0 x0166 y0166 z0166 M  
 C-CT--0.319 0 x0167 y0167 z0167 M  
 H-H-0.219 0 x0168 y0168 z0168 M  
 H-H1-0.076 0 x0169 y0169 z0169 M  
 H-HC--0.030 0 x0170 y0170 z0170 M  
 H-HC-0.079 0 x0171 y0171 z0171 M  
 H-HC-0.079 0 x0172 y0172 z0172 M  
 H-HC-0.079 0 x0173 y0173 z0173 M  
 H-HC-0.079 0 x0174 y0174 z0174 M  
 H-HC-0.079 0 x0175 y0175 z0175 M  
 H-HC-0.079 0 x0176 y0176 z0176 M  
 N-N--0.416 0 x0177 y0177 z0177 M  
 C-CT--0.059 0 x0178 y0178 z0178 M  
 C-C-0.337 0 x0179 y0179 z0179 M  
 O-O--0.468 0 x0180 y0180 z0180 M  
 C-CT-0.130 0 x0181 y0181 z0181 M  
 C-CT--0.043 0 x0182 y0182 z0182 M  
 C-CT--0.320 0 x0183 y0183 z0183 M  
 C-CT--0.066 0 x0184 y0184 z0184 M  
 H-H-0.272 0 x0185 y0185 z0185 M  
 H-H1-0.056 0 x0186 y0186 z0186 M  
 H-HC-0.019 0 x0187 y0187 z0187 M  
 H-HC-0.024 0 x0188 y0188 z0188 M  
 H-HC-0.024 0 x0189 y0189 z0189 M  
 H-HC-0.088 0 x0190 y0190 z0190 M  
 H-HC-0.088 0 x0191 y0191 z0191 M  
 H-HC-0.086 0 x0192 y0192 z0192 M  
 H-HC-0.019 0 x0193 y0193 z0193 M  
 H-HC-0.019 0 x0194 y0194 z0194 M  
 H-HC-0.019 0 x0195 y0195 z0195 M  
 N-N--0.416 0 x0196 y0196 z0196 M  
 C-CT--0.087 0 x0197 y0197 z0197 M  
 C-C-0.597 0 x0198 y0198 z0198 M  
 O-O--0.568 0 x0199 y0199 z0199 M  
 C-CT-0.152 0 x0200 y0200 z0200 M  
 C-CT--0.148 0 x0201 y0201 z0201 M  
 C-CT--0.059 0 x0202 y0202 z0202 M  
 H-H-0.267 0 x0203 y0203 z0203 M  
 H-H1-0.097 0 x0204 y0204 z0204 M  
 H-HC--0.011 0 x0205 y0205 z0205 M  
 H-HC-0.072 0 x0206 y0206 z0206 M  
 H-HC-0.031 0 x0207 y0207 z0207 M  
 H-HC-0.015 0 x0208 y0208 z0208 M  
 H-HC-0.008 0 x0209 y0209 z0209 M  
 H-HC-0.017 0 x0210 y0210 z0210 M  
 H-HC-0.013 0 x0211 y0211 z0211 M  
 C-CT--0.087 0 x0212 y0212 z0212 M  
 O-O--0.568 0 x0213 y0213 z0213 M  
 H-N--0.416 0 x0214 y0214 z0214 M  
 C-C-0.597 0 x0215 y0215 z0215 M  
 H-CT-0.298 0 x0216 y0216 z0216 M  
 H-H1-0.097 0 x0217 y0217 z0217 M  
 C-CT-0.034 0 x0218 y0218 z0218 M  
 H-H-0.272 0 x0219 y0219 z0219 M  
 N-N--0.416 0 x0220 y0220 z0220 M  
 H-C-0.412 0 x0221 y0221 z0221 M  
 H-CT--0.182 0 x0222 y0222 z0222 M  
 H-H1-0.082 0 x0223 y0223 z0223 M  
 C-CT--0.022 0 x0224 y0224 z0224 M  
 O-O--0.133 0 x0225 y0225 z0225 M  
 H-N--0.121 0 x0226 y0226 z0226 M  
 C-C-0.203 0 x0227 y0227 z0227 M  
 H-CT-0.034 0 x0228 y0228 z0228 M  
 H-H1-0.037 0 x0229 y0229 z0229 M  
 C-CT-0.016 0 x0230 y0230 z0230 M  
 H-H-0.241 0 x0231 y0231 z0231 M  
 N-N--0.368 0 x0232 y0232 z0232 M  
 H-C-0.241 0 x0233 y0233 z0233 M  
 H-CT--0.034 0 x0234 y0234 z0234 M  
 H-H1-0.050 0 x0235 y0235 z0235 M  
 C-CT--0.023 0 x0236 y0236 z0236 M  
 O-O--0.382 0 x0237 y0237 z0237 M  
 H-N3-0.074 0 x0238 y0238 z0238 M  
 C-C-0.597 0 x0239 y0239 z0239 M  
 H-CT-0.025 0 x0240 y0240 z0240 M  
 H-H1-0.088 0 x0241 y0241 z0241 M  
 C-CT--0.087 0 x0242 y0242 z0242    
 H-H-0.272 0 x0243 y0243 z0243    
 N-N--0.416 0 x0244 y0244 z0244    
 H-C-0.597 0 x0245 y0245 z0245    
 H-CT-0.298 0 x0246 y0246 z0246    
 H-H1-0.097 0 x0247 y0247 z0247    
 C-CT--0.016 0 x0248 y0248 z0248 M  
 O-O--0.286 0 x0249 y0249 z0249 M  
 H-N3-0.175 0 x0250 y0250 z0250 M  
 C-C-0.374 0 x0251 y0251 z0251 M  
 H-CT-0.009 0 x0252 y0252 z0252 M  
 H-H1-0.063 0 x0253 y0253 z0253 M  
 C-CT--0.009 0 x0254 y0254 z0254 M  
 H-H-0.211 0 x0255 y0255 z0255 M  
 N-N--0.242 0 x0256 y0256 z0256 M  
 H-C-0.355 0 x0257 y0257 z0257 M  
 H-H1-0.011 0 x0258 y0258 z0258 M  
 H-H1-0.035 0 x0259 y0259 z0259 M  
 C-CT--0.003 0 x0260 y0260 z0260 M  
 O-O--0.114 0 x0261 y0261 z0261 M  
 H-N3-0.024 0 x0262 y0262 z0262 M  
 C-C-0.188 0 x0263 y0263 z0263 M  
 H-CT-0.003 0 x0264 y0264 z0264 M  
 H-H1-0.016 0 x0265 y0265 z0265 M  
 C-CT--0.025 0 x0266 y0266 z0266 M  
 H-H-0.272 0 x0267 y0267 z0267 M  
 N-N--0.416 0 x0268 y0268 z0268 M  
 H-C-0.597 0 x0269 y0269 z0269 M  
 H-H1-0.061 0 x0270 y0270 z0270 M  
 H-H1-0.070 0 x0271 y0271 z0271 M  
 C-CT--0.087 0 x0272 y0272 z0272 M  
 O-O--0.568 0 x0273 y0273 z0273 M  
 H-N--0.416 0 x0274 y0274 z0274 M  
 C-C-0.597 0 x0275 y0275 z0275 M  
 H-CT-0.298 0 x0276 y0276 z0276 M  
 H-H1-0.097 0 x0277 y0277 z0277 M  
 C-CT--0.087 0 x0278 y0278 z0278 M  
 H-H-0.272 0 x0279 y0279 z0279 M  
 N-N--0.416 0 x0280 y0280 z0280 M  
 H-C-0.597 0 x0281 y0281 z0281 M  
 H-CT-0.298 0 x0282 y0282 z0282 M  
 H-H1-0.097 0 x0283 y0283 z0283 M  
 C-CT--0.005 0 x0284 y0284 z0284 M  
 O-O--0.119 0 x0285 y0285 z0285 M  
 H-N--0.106 0 x0286 y0286 z0286 M  
 C-C-0.180 0 x0287 y0287 z0287 M  
 H-H1-0.007 0 x0288 y0288 z0288 M  
 H-H1-0.022 0 x0289 y0289 z0289 M  
 C-CT-0.003 0 x0290 y0290 z0290 M  
 H-H-0.047 0 x0291 y0291 z0291 M  
 N-N--0.073 0 x0292 y0292 z0292 M  
 H-C-0.032 0 x0293 y0293 z0293 M  
 H-CT--0.014 0 x0294 y0294 z0294 M  
 H-H1-0.010 0 x0295 y0295 z0295 M  
 Variables:
x0008 =    -1.4346
y0008 =     3.3364
z0008 =    15.9049
x0009 =     1.1894
y0009 =     2.7905
z0009 =    17.2981
x0010 =     0.1297
y0010 =     3.4623
z0010 =    19.4660
x0011 =     1.3700
y0011 =     5.2539
z0011 =    18.6621
x0012 =    -0.2651
y0012 =     5.9985
z0012 =    18.6621
x0013 =     0.4278
y0013 =     5.3524
z0013 =    17.1355
x0014 =    -2.1476
y0014 =     3.7332
z0014 =    19.5693
x0015 =    -2.2585
y0015 =     3.4323
z0015 =    17.8015
x0016 =    -2.0327
y0016 =     5.1056
z0016 =    18.4154
 Constants:
x0001 =    -0.4610
y0001 =     3.4040
z0001 =    16.1650
x0002 =     0.1000
y0002 =     2.9010
z0002 =    17.4030
x0003 =    -0.5380
y0003 =     1.5840
z0003 =    17.7320
x0004 =    -1.7120
y0004 =     1.4560
z0004 =    17.5310
x0005 =    -0.2500
y0005 =     3.8640
z0005 =    18.5150
x0006 =     0.3620
y0006 =     5.2080
z0006 =    18.2240
x0007 =    -1.7830
y0007 =     4.0470
z0007 =    18.5800
x0017 =    -0.4240
y0017 =     8.4720
z0017 =    11.1930
x0018 =     0.1590
y0018 =     7.9090
z0018 =    12.4230
x0019 =    -0.4980
y0019 =     6.6220
z0019 =    12.8910
x0020 =    -1.6670
y0020 =     6.3860
z0020 =    12.6440
x0021 =    -1.4146
y0021 =     8.5333
z0021 =    11.0055
x0022 =     0.0622
y0022 =     8.6566
z0022 =    13.2240
x0023 =     1.2239
y0023 =     7.7050
z0023 =    12.2378
x0024 =     0.2370
y0024 =     5.7780
z0024 =    13.5850
x0025 =    -0.3900
y0025 =     4.5590
z0025 =    14.0790
x0026 =     0.3180
y0026 =     4.0000
z0026 =    15.2850
x0027 =     1.5440
y0027 =     4.0730
z0027 =    15.4150
x0028 =     1.2012
y0028 =     6.0435
z0028 =    13.7266
x0029 =    -0.3749
y0029 =     3.8044
z0029 =    13.2788
x0030 =    -1.4315
y0030 =     4.7831
z0030 =    14.3528
x0031 =     0.2190
y0031 =     0.6320
z0031 =    18.2550
x0032 =    -0.3300
y0032 =    -0.6090
z0032 =    18.8190
x0033 =     0.2210
y0033 =    -0.8300
z0033 =    20.2020
x0034 =     1.3840
y0034 =    -0.4990
z0034 =    20.4980
x0035 =     0.0720
y0035 =    -1.8740
z0035 =    17.9890
x0036 =     1.4760
y0036 =    -2.1950
z0036 =    18.2740
x0037 =    -0.7130
y0037 =    -3.0980
z0037 =    18.3670
x0038 =     1.2047
y0038 =     0.8501
z0038 =    18.2241
x0039 =    -1.4233
y0039 =    -0.4878
z0039 =    18.8136
x0040 =    -0.1195
y0040 =    -1.6304
z0040 =    16.9335
x0041 =     1.5626
y0041 =    -2.5787
z0041 =    19.3013
x0042 =     2.0879
y0042 =    -1.2870
z0042 =    18.1684
x0043 =     1.8279
y0043 =    -2.9594
z0043 =    17.5656
x0044 =    -1.7886
y0044 =    -2.8830
z0044 =    18.2848
x0045 =    -0.4733
y0045 =    -3.3814
z0045 =    19.4025
x0046 =    -0.4510
y0046 =    -3.9245
z0046 =    17.6901
x0047 =    -0.6170
y0047 =    -1.4240
z0047 =    21.0540
x0048 =    -0.2040
y0048 =    -1.9670
z0048 =    22.3760
x0049 =    -0.7150
y0049 =    -3.3900
z0049 =    22.4610
x0050 =    -1.8590
y0050 =    -3.6610
z0050 =    22.1310
x0051 =    -0.8210
y0051 =    -1.2060
z0051 =    23.5590
x0052 =    -0.5010
y0052 =     0.2810
z0052 =    23.5050
x0053 =    -0.3120
y0053 =    -1.8050
z0053 =    24.8930
x0054 =    -0.9410
y0054 =     1.0210
z0054 =    24.7230
x0055 =    -1.5618
y0055 =    -1.4625
z0055 =    20.6990
x0056 =     0.8906
y0056 =    -1.8817
z0056 =    22.4431
x0057 =    -1.9136
y0057 =    -1.3160
z0057 =    23.4941
x0058 =    -1.0066
y0058 =     0.7158
z0058 =    22.6302
x0059 =     0.5872
y0059 =     0.4007
z0059 =    23.3973
x0060 =     0.4863
y0060 =    -1.1669
z0060 =    25.2999
x0061 =     0.0818
y0061 =    -2.8163
z0061 =    24.7136
x0062 =    -1.1430
y0062 =    -1.8567
z0062 =    25.6118
x0063 =    -0.4679
y0063 =     0.5777
z0063 =    25.6116
x0064 =    -2.0349
y0064 =     0.9535
z0064 =    24.8172
x0065 =    -0.6443
y0065 =     2.0767
z0065 =    24.6371
x0066 =     1.8990
y0066 =    -1.3840
z0066 =    11.1140
x0067 =     2.2910
y0067 =    -2.0300
z0067 =    12.4000
x0068 =     1.7300
y0068 =    -3.4450
z0068 =    12.4590
x0069 =     0.5170
y0069 =    -3.6570
z0069 =    12.4520
x0070 =     1.7560
y0070 =    -1.3000
z0070 =    13.6290
x0071 =     2.1440
y0071 =     0.1860
z0071 =    13.6130
x0072 =     2.2700
y0072 =    -2.0380
z0072 =    14.8890
x0073 =     1.5720
y0073 =     1.0030
z0073 =    14.7320
x0074 =     0.9402
y0074 =    -1.2906
z0074 =    10.8107
x0075 =     3.3906
y0075 =    -2.0104
z0075 =    12.4204
x0076 =     0.6561
y0076 =    -1.3147
z0076 =    13.6310
x0077 =     1.7971
y0077 =     0.6187
z0077 =    12.6630
x0078 =     3.2406
y0078 =     0.2527
z0078 =    13.6682
x0079 =     2.3993
y0079 =    -1.3162
z0079 =    15.7089
x0080 =     3.2350
y0080 =    -2.5161
z0080 =    14.6651
x0081 =     1.5404
y0081 =    -2.8054
z0081 =    15.1870
x0082 =     1.2263
y0082 =     0.3342
z0082 =    15.5340
x0083 =     0.7243
y0083 =     1.5952
z0083 =    14.3568
x0084 =     2.3462
y0084 =     1.6783
z0084 =    15.1252
x0085 =    -4.5770
y0085 =    -2.6960
z0085 =    22.2940
x0086 =    -5.1330
y0086 =    -1.4370
z0086 =    21.9330
x0087 =    -4.3810
y0087 =    -0.8860
z0087 =    20.6730
x0088 =    -3.1570
y0088 =    -0.8090
z0088 =    20.6610
x0089 =    -4.9660
y0089 =    -0.4400
z0089 =    23.1630
x0090 =    -5.4420
y0090 =     1.0120
z0090 =    22.8590
x0091 =    -5.5690
y0091 =    -0.9980
z0091 =    24.4150
x0092 =    -3.5902
y0092 =    -2.8727
z0092 =    22.4173
x0093 =    -6.2009
y0093 =    -1.5387
z0093 =    21.6895
x0094 =    -3.8837
y0094 =    -0.3510
z0094 =    23.3382
x0095 =    -6.2357
y0095 =     1.2929
z0095 =    23.5669
x0096 =    -5.8317
y0096 =     1.0620
z0096 =    21.8316
x0097 =    -4.5941
y0097 =     1.7049
z0097 =    22.9642
x0098 =    -4.8367
y0098 =    -1.6484
z0098 =    24.9158
x0099 =    -6.4656
y0099 =    -1.5823
z0099 =    24.1607
x0100 =    -5.8470
y0100 =    -0.1726
z0100 =    25.0869
x0101 =    -5.1300
y0101 =    -0.5090
z0101 =    19.6350
x0102 =    -4.5770
y0102 =     0.0730
z0102 =    18.4260
x0103 =    -5.2850
y0103 =     1.3320
z0103 =    17.9870
x0104 =    -6.4970
y0104 =     1.4790
z0104 =    18.1610
x0105 =    -6.1194
y0105 =    -0.6648
z0105 =    19.7649
x0106 =    -4.6513
y0106 =    -0.6685
z0106 =    17.6169
x0107 =    -3.5191
y0107 =     0.3136
z0107 =    18.6076
x0108 =    -4.5380
y0108 =     2.2570
z0108 =    17.3990
x0109 =    -5.1550
y0109 =     3.4470
z0109 =    16.8150
x0110 =    -4.4330
y0110 =     3.9040
z0110 =    15.5600
x0111 =    -3.1960
y0111 =     3.7410
z0111 =    15.4470
x0112 =    -3.5465
y0112 =     2.0644
z0112 =    17.3961
x0113 =    -5.1291
y0113 =     4.2600
z0113 =    17.5555
x0114 =    -6.1997
y0114 =     3.2162
z0114 =    16.5595
x0115 =    -2.1560
y0115 =    -2.5460
z0115 =    12.1250
x0116 =    -2.6900
y0116 =    -1.2810
z0116 =    11.7490
x0117 =    -1.9340
y0117 =    -0.8440
z0117 =    10.4470
x0118 =    -0.7010
y0118 =    -1.0090
z0118 =    10.3210
x0119 =    -2.4690
y0119 =    -0.2920
z0119 =    12.9080
x0120 =    -2.9720
y0120 =     1.1210
z0120 =    12.5760
x0121 =    -3.1000
y0121 =    -0.8160
z0121 =    14.1700
x0122 =    -1.1743
y0122 =    -2.7065
z0122 =    12.3001
x0123 =    -3.7709
y0123 =    -1.3162
z0123 =    11.5477
x0124 =    -1.3840
y0124 =    -0.2058
z0124 =    13.0670
x0125 =    -3.7776
y0125 =     1.3960
z0125 =    13.2727
x0126 =    -3.3555
y0126 =     1.1404
z0126 =    11.5452
x0127 =    -2.1426
y0127 =     1.8370
z0127 =    12.6727
x0128 =    -2.8514
y0128 =    -0.1484
z0128 =    15.0082
x0129 =    -2.7171
y0129 =    -1.8258
z0129 =    14.3788
x0130 =    -4.1920
y0130 =    -0.8562
z0130 =    14.0440
x0131 =    -2.1710
y0131 =    10.8300
z0131 =    18.7130
x0132 =    -2.6450
y0132 =     9.4710
z0132 =    18.8160
x0133 =    -1.9450
y0133 =     8.9880
z0133 =    20.0900
x0134 =    -0.7300
y0134 =     9.1980
z0134 =    20.2400
x0135 =    -2.1980
y0135 =     8.5510
z0135 =    17.5770
x0136 =    -2.8030
y0136 =     7.1290
z0136 =    17.6500
x0137 =    -2.5190
y0137 =     9.1490
z0137 =    16.2230
x0138 =    -1.1901
y0138 =    11.0705
z0138 =    18.7199
x0139 =    -3.7435
y0139 =     9.4154
z0139 =    18.8272
x0140 =    -1.1037
y0140 =     8.4910
z0140 =    17.6709
x0141 =    -2.6017
y0141 =     6.6959
z0141 =    18.6409
x0142 =    -2.3473
y0142 =     6.4969
z0142 =    16.8736
x0143 =    -3.8895
y0143 =     7.1846
z0143 =    17.4874
x0144 =    -1.6204
y0144 =     9.6348
z0144 =    15.8150
x0145 =    -3.3213
y0145 =     9.8935
z0145 =    16.3325
x0146 =    -2.8481
y0146 =     8.3523
z0146 =    15.5396
x0147 =     4.3400
y0147 =    -0.5090
z0147 =    19.6350
x0148 =     4.8930
y0148 =     0.0730
z0148 =    18.4260
x0149 =     4.1850
y0149 =     1.3320
z0149 =    17.9870
x0150 =     2.9730
y0150 =     1.4790
z0150 =    18.1610
x0151 =     3.3506
y0151 =    -0.6648
z0151 =    19.7649
x0152 =     4.8187
y0152 =    -0.6685
z0152 =    17.6169
x0153 =     5.9509
y0153 =     0.3136
z0153 =    18.6076
x0154 =     4.9320
y0154 =     2.2570
z0154 =    17.3990
x0155 =     4.3150
y0155 =     3.4470
z0155 =    16.8150
x0156 =     5.0370
y0156 =     3.9040
z0156 =    15.5600
x0157 =     6.2740
y0157 =     3.7410
z0157 =    15.4470
x0158 =     5.9235
y0158 =     2.0644
z0158 =    17.3961
x0159 =     4.3409
y0159 =     4.2600
z0159 =    17.5555
x0160 =     3.2703
y0160 =     3.2162
z0160 =    16.5595
x0161 =     1.7630
y0161 =     4.6290
z0161 =    24.6070
x0162 =     2.4350
y0162 =     5.2040
z0162 =    23.4390
x0163 =     1.7620
y0163 =     6.4460
z0163 =    22.9300
x0164 =     0.5430
y0164 =     6.4960
z0164 =    22.8520
x0165 =     2.4770
y0165 =     4.1600
z0165 =    22.2740
x0166 =     1.1590
y0166 =     3.5340
z0166 =    22.0970
x0167 =     3.0260
y0167 =     4.8060
z0167 =    20.8840
x0168 =     0.7696
y0168 =     4.6750
z0168 =    24.7832
x0169 =     3.4496
y0169 =     5.4720
z0169 =    23.7690
x0170 =     3.1965
y0170 =     3.3771
z0170 =    22.5559
x0171 =     0.6423
y0171 =     3.4880
z0171 =    23.0670
x0172 =     1.2874
y0172 =     2.5159
z0172 =    21.7007
x0173 =     0.5623
y0173 =     4.1303
z0173 =    21.3910
x0174 =     3.8515
y0174 =     5.4992
z0174 =    21.1033
x0175 =     2.2106
y0175 =     5.3517
z0175 =    20.3867
x0176 =     3.3863
y0176 =     4.0030
z0176 =    20.2242
x0177 =     1.8020
y0177 =     9.5740
z0177 =    19.5860
x0178 =     2.2160
y0178 =    10.4740
z0178 =    18.4940
x0179 =     1.2950
y0179 =    11.6670
z0179 =    18.5970
x0180 =     0.0810
y0180 =    11.4480
z0180 =    18.6050
x0181 =     1.8840
y0181 =     9.9090
z0181 =    17.0970
x0182 =     2.3830
y0182 =     8.4690
z0182 =    16.9370
x0183 =     2.4530
y0183 =    10.8400
z0183 =    16.0150
x0184 =     1.5310
y0184 =     7.6140
z0184 =    16.0460
x0185 =     0.8350
y0185 =     9.3531
z0185 =    19.7765
x0186 =     3.2966
y0186 =    10.6566
z0186 =    18.5883
x0187 =     0.7908
y0187 =     9.8702
z0187 =    16.9809
x0188 =     2.4147
y0188 =     8.0032
z0188 =    17.9330
x0189 =     3.3981
y0189 =     8.5016
z0189 =    16.5144
x0190 =     1.8872
y0190 =    11.7833
z0190 =    16.0091
x0191 =     2.3672
y0191 =    10.3542
z0191 =    15.0318
x0192 =     3.5114
y0192 =    11.0481
z0192 =    16.2306
x0193 =     0.7673
y0193 =     7.1027
z0193 =    16.6504
x0194 =     2.1622
y0194 =     6.8666
z0194 =    15.5431
x0195 =     1.0397
y0195 =     8.2467
z0195 =    15.2921
x0196 =     2.5550
y0196 =     7.4390
z0196 =    22.5550
x0197 =     2.0380
y0197 =     8.5730
z0197 =    21.7950
x0198 =     2.6620
y0198 =     8.9650
z0198 =    20.4270
x0199 =     3.8650
y0199 =     8.7310
z0199 =    20.1490
x0200 =     1.8890
y0200 =     9.8710
z0200 =    22.6810
x0201 =     0.5430
y0201 =     9.8730
z0201 =    23.2720
x0202 =     2.9270
y0202 =     9.9360
z0202 =    23.7690
x0203 =     3.5189
y0203 =     7.3353
z0203 =    22.8383
x0204 =     1.0715
y0204 =     8.1380
z0204 =    21.5005
x0205 =     2.0383
y0205 =    10.7534
z0205 =    22.0415
x0206 =     0.0637
y0206 =     8.8990
z0206 =    23.0944
x0207 =    -0.0591
y0207 =    10.6685
z0207 =    22.8087
x0208 =     0.6187
y0208 =    10.0530
z0208 =    24.3545
x0209 =     2.9857
y0209 =    10.9632
z0209 =    24.1581
x0210 =     3.9049
y0210 =     9.6418
z0210 =    23.3601
x0211 =     2.6479
y0211 =     9.2505
z0211 =    24.5827
x0212 =    -0.3470
y0212 =     9.4920
z0212 =     8.9950
x0213 =     1.5840
y0213 =     8.9700
z0213 =    10.2830
x0214 =     0.0644
y0214 =    10.5121
z0214 =     9.0043
x0215 =     0.3510
y0215 =     8.9510
z0215 =    10.2200
x0216 =    -0.0165
y0216 =     8.9314
z0216 =     8.1081
x0217 =    -1.4454
y0217 =     9.4345
z0217 =     8.9821
x0218 =    -0.3740
y0218 =    -5.6570
z0218 =    22.8580
x0219 =     0.9953
y0219 =    -4.0038
z0219 =    23.2395
x0220 =     0.0940
y0220 =    -4.3080
z0220 =    22.9000
x0221 =     0.2920
y0221 =    -6.3564
z0221 =    23.3846
x0222 =    -0.4785
y0222 =    -5.9849
z0222 =    21.8132
x0223 =    -1.3471
y0223 =    -5.6691
z0223 =    23.3708
x0224 =     2.1700
y0224 =    -0.5480
z0224 =     8.8500
x0225 =     3.9570
y0225 =    -0.6360
z0225 =    10.4760
x0226 =     2.5059
y0226 =     0.4220
z0226 =     8.4548
x0227 =     2.7590
y0227 =    -0.8600
z0227 =    10.2240
x0228 =     2.4558
y0228 =    -1.3662
z0228 =     8.1726
x0229 =     1.0761
y0229 =    -0.4718
z0229 =     8.9368
x0230 =     2.1020
y0230 =    -5.7450
z0230 =    12.2280
x0231 =     3.5527
y0231 =    -4.1767
z0231 =    12.6414
x0232 =     2.5890
y0232 =    -4.4230
z0232 =    12.4660
x0233 =     2.4996
y0233 =    -6.4746
z0233 =    12.9489
x0234 =     2.3433
y0234 =    -6.0490
z0234 =    11.1987
x0235 =     1.0104
y0235 =    -5.7280
z0235 =    12.3628
x0236 =    -4.6760
y0236 =    -4.9640
z0236 =    23.0200
x0237 =    -6.5700
y0237 =    -3.6770
z0237 =    22.2760
x0238 =    -5.2259
y0238 =    -5.2255
z0238 =    23.9361
x0239 =    -5.3500
y0239 =    -3.7170
z0239 =    22.4890
x0240 =    -4.7523
y0240 =    -5.8010
z0240 =    22.3104
x0241 =    -3.6053
y0241 =    -4.7795
z0241 =    23.1922
x0242 =    -4.6080
y0242 =     4.9660
z0242 =    13.3390
x0243 =    -6.1656
y0243 =     4.5575
z0243 =    14.8868
x0244 =    -5.1920
y0244 =     4.4850
z0244 =    14.6280
x0245 =    -5.0364
y0245 =     5.9242
z0245 =    13.0099
x0246 =    -4.7251
y0246 =     4.1957
z0246 =    12.5625
x0247 =    -3.5370
y0247 =     5.1484
z0247 =    13.5113
x0248 =    -2.1160
y0248 =    -4.8860
z0248 =    12.3220
x0249 =    -4.1290
y0249 =    -3.5680
z0249 =    12.3980
x0250 =    -2.2211
y0250 =    -5.3286
z0250 =    13.3235
x0251 =    -2.9000
y0251 =    -3.6020
z0251 =    12.2730
x0252 =    -2.4235
y0252 =    -5.5679
z0252 =    11.5155
x0253 =    -1.0485
y0253 =    -4.6843
z0253 =    12.1495
x0254 =    -2.1860
y0254 =     0.2550
z0254 =     8.2550
x0255 =    -3.6923
y0255 =    -0.4051
z0255 =     9.6587
x0256 =    -2.7020
y0256 =    -0.3380
z0256 =     9.4720
x0257 =    -2.6772
y0257 =     1.2055
z0257 =     7.9996
x0258 =    -2.3401
y0258 =    -0.4529
z0258 =     7.4272
x0259 =    -1.1094
y0259 =     0.4402
z0259 =     8.3842
x0260 =    -2.3230
y0260 =    13.2420
z0260 =    18.7140
x0261 =    -4.1940
y0261 =    11.7740
z0261 =    18.3890
x0262 =    -2.6474
y0262 =    13.8784
z0262 =    17.8775
x0263 =    -2.9820
y0263 =    11.8740
z0263 =    18.6010
x0264 =    -2.5597
y0264 =    13.7238
z0264 =    19.6741
x0265 =    -1.2324
y0265 =    13.1053
z0265 =    18.6698
x0266 =    -2.0970
y0266 =     7.8530
z0266 =    22.2630
x0267 =    -3.6438
y0267 =     8.1738
z0267 =    20.7186
x0268 =    -2.6820
y0268 =     8.3200
z0268 =    20.9900
x0269 =    -2.6582
y0269 =     7.0393
z0269 =    22.7456
x0270 =    -2.0587
y0270 =     8.7048
z0270 =    22.9579
x0271 =    -1.0759
y0271 =     7.4962
z0271 =    22.0629
x0272 =     4.3370
y0272 =    -1.4370
z0272 =    21.9330
x0273 =     6.3130
y0273 =    -0.8090
z0273 =    20.6610
x0274 =     4.7668
y0274 =    -2.4103
z0274 =    22.2121
x0275 =     5.0890
y0275 =    -0.8860
z0275 =    20.6730
x0276 =     4.4524
y0276 =    -0.7482
z0276 =    22.7828
x0277 =     3.2691
y0277 =    -1.5387
z0277 =    21.6895
x0278 =     4.8620
y0278 =     4.9660
z0278 =    13.3390
x0279 =     3.3044
y0279 =     4.5575
z0279 =    14.8868
x0280 =     4.2780
y0280 =     4.4850
z0280 =    14.6280
x0281 =     4.4336
y0281 =     5.9242
z0281 =    13.0099
x0282 =     4.7449
y0282 =     4.1957
z0282 =    12.5625
x0283 =     5.9330
y0283 =     5.1484
z0283 =    13.5113
x0284 =     1.7180
y0284 =     3.5230
z0284 =    26.7350
x0285 =     3.7030
y0285 =     3.7900
z0285 =    25.4440
x0286 =     2.2337
y0286 =     2.6635
z0286 =    27.1881
x0287 =     2.4720
y0287 =     3.9900
z0287 =    25.5300
x0288 =     1.6589
y0288 =     4.3410
z0288 =    27.4681
x0289 =     0.7026
y0289 =     3.2233
z0289 =    26.4365
x0290 =     2.9690
y0290 =    13.3180
z0290 =    19.3950
x0291 =     1.1210
y0291 =    13.5845
z0291 =    18.3538
x0292 =     1.7510
y0292 =    12.9040
z0292 =    18.7540
x0293 =     2.5668
y0293 =    13.8133
z0293 =    20.2910
x0294 =     3.4732
y0294 =    14.0969
z0294 =    18.8042
x0295 =     3.6821
y0295 =    12.4957
z0295 =    19.5542
 
HrmStr1   * * 331. 1.09
HrmBnd1   * * * 20. 90.
AmbTrs    * * * * 0 0 0 0 0.0 0.0 1.15 0.0 3.0
VDW       MG 1.09 0.25
VDW       NT 1.8240 0.17

   1-0295 0
6-31G*  
****
 
END

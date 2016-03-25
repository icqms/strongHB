g09 << +++ > VAL5-B.log
%mem=1900MB
%Nprocshared=4
%Chk=/short/d63/reimers/VAL5-B
#P ONIOM(cam-b3lyp:AMBER=hardfirst)=EmbedCharge iop(4/28=5)
   opt(z-matrix,maxcycle=7,nrscale) nosym scf(conver=6) 

  lysozyme partial opt VAL5-B

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
 N-N--0.414 0 x0017 y0017 z0017 M  
 C-CT--0.087 0 x0018 y0018 z0018 M  
 C-C-0.597 0 x0019 y0019 z0019 M  
 O-O--0.568 0 x0020 y0020 z0020 M  
 C-CT-0.298 0 x0021 y0021 z0021 M  
 C-CT--0.319 0 x0022 y0022 z0022 M  
 C-CT--0.319 0 x0023 y0023 z0023 M  
 H-H-0.269 0 x0024 y0024 z0024 M  
 H-H1-0.097 0 x0025 y0025 z0025 M  
 H-HC--0.030 0 x0026 y0026 z0026 M  
 H-HC-0.079 0 x0027 y0027 z0027 M  
 H-HC-0.079 0 x0028 y0028 z0028 M  
 H-HC-0.079 0 x0029 y0029 z0029 M  
 H-HC-0.079 0 x0030 y0030 z0030 M  
 H-HC-0.077 0 x0031 y0031 z0031 M  
 H-HC-0.079 0 x0032 y0032 z0032 M  
 N-N--0.416 0 x0033 y0033 z0033 M H-H1--0.144 
 C-CT--0.025 0 x0034 y0034 z0034    
 C-C-0.597 0 x0035 y0035 z0035    
 O-O--0.568 0 x0036 y0036 z0036    
 H-H-0.272 0 x0037 y0037 z0037 M  
 H-H1-0.070 0 x0038 y0038 z0038    
 H-H1-0.070 0 x0039 y0039 z0039    
 N-N--0.416 0 x0040 y0040 z0040    
 C-CT--0.025 0 x0041 y0041 z0041    
 C-C-0.597 0 x0042 y0042 z0042    
 O-O--0.568 0 x0043 y0043 z0043    
 H-H-0.272 0 x0044 y0044 z0044    
 H-H1-0.070 0 x0045 y0045 z0045    
 H-H1-0.070 0 x0046 y0046 z0046    
 N-N--0.416 0 x0047 y0047 z0047    
 C-CT--0.087 0 x0048 y0048 z0048    
 C-C-0.597 0 x0049 y0049 z0049    
 O-O--0.568 0 x0050 y0050 z0050    
 C-CT-0.298 0 x0051 y0051 z0051    
 C-CT--0.319 0 x0052 y0052 z0052    
 C-CT--0.319 0 x0053 y0053 z0053    
 H-H-0.272 0 x0054 y0054 z0054    
 H-H1-0.097 0 x0055 y0055 z0055    
 H-HC--0.030 0 x0056 y0056 z0056    
 H-HC-0.079 0 x0057 y0057 z0057    
 H-HC-0.079 0 x0058 y0058 z0058    
 H-HC-0.079 0 x0059 y0059 z0059    
 H-HC-0.079 0 x0060 y0060 z0060    
 H-HC-0.079 0 x0061 y0061 z0061    
 H-HC-0.079 0 x0062 y0062 z0062    
 N-N--0.416 0 x0063 y0063 z0063    
 C-CT--0.060 0 x0064 y0064 z0064    
 C-C-0.597 0 x0065 y0065 z0065 M H-H1-0.058  
 O-O--0.568 0 x0066 y0066 z0066 M  
 C-CT-0.130 0 x0067 y0067 z0067 M H-H1-0.058  
 C-CT--0.043 0 x0068 y0068 z0068 M  
 C-CT--0.320 0 x0069 y0069 z0069 M  
 C-CT--0.066 0 x0070 y0070 z0070 M  
 H-H-0.272 0 x0071 y0071 z0071    
 H-H1-0.087 0 x0072 y0072 z0072    
 H-HC-0.019 0 x0073 y0073 z0073 M  
 H-HC-0.024 0 x0074 y0074 z0074 M  
 H-HC-0.024 0 x0075 y0075 z0075 M  
 H-HC-0.088 0 x0076 y0076 z0076 M  
 H-HC-0.088 0 x0077 y0077 z0077 M  
 H-HC-0.088 0 x0078 y0078 z0078 M  
 H-HC-0.019 0 x0079 y0079 z0079 M  
 H-HC-0.019 0 x0080 y0080 z0080 M  
 H-HC-0.019 0 x0081 y0081 z0081 M  
 N-N--0.416 0 x0082 y0082 z0082 M  
 C-CT--0.087 0 x0083 y0083 z0083 M  
 C-C-0.597 0 x0084 y0084 z0084 M  
 O-O--0.568 0 x0085 y0085 z0085 M  
 C-CT-0.298 0 x0086 y0086 z0086 M  
 C-CT--0.319 0 x0087 y0087 z0087 M  
 C-CT--0.297 0 x0088 y0088 z0088 M  
 H-H-0.272 0 x0089 y0089 z0089 M  
 H-H1-0.097 0 x0090 y0090 z0090 M  
 H-HC--0.030 0 x0091 y0091 z0091 M  
 H-HC-0.079 0 x0092 y0092 z0092 M  
 H-HC-0.079 0 x0093 y0093 z0093 M  
 H-HC-0.079 0 x0094 y0094 z0094 M  
 H-HC-0.050 0 x0095 y0095 z0095 M  
 H-HC-0.078 0 x0096 y0096 z0096 M  
 H-HC-0.036 0 x0097 y0097 z0097 M  
 N-N--0.416 0 x0098 y0098 z0098 M H-H1--0.144 
 C-CT--0.025 0 x0099 y0099 z0099    
 C-C-0.597 0 x0100 y0100 z0100    
 O-O--0.568 0 x0101 y0101 z0101    
 H-H-0.272 0 x0102 y0102 z0102 M  
 H-H1-0.070 0 x0103 y0103 z0103    
 H-H1-0.070 0 x0104 y0104 z0104    
 N-N--0.416 0 x0105 y0105 z0105    
 C-CT--0.025 0 x0106 y0106 z0106    
 C-C-0.597 0 x0107 y0107 z0107    
 O-O--0.568 0 x0108 y0108 z0108    
 H-H-0.272 0 x0109 y0109 z0109    
 H-H1-0.070 0 x0110 y0110 z0110    
 H-H1-0.070 0 x0111 y0111 z0111    
 N-N--0.213 0 x0112 y0112 z0112 M  
 C-CT--0.060 0 x0113 y0113 z0113 M  
 C-C-0.597 0 x0114 y0114 z0114 M  
 O-O--0.469 0 x0115 y0115 z0115 M  
 C-CT-0.130 0 x0116 y0116 z0116 M  
 C-CT--0.043 0 x0117 y0117 z0117 M  
 C-CT--0.320 0 x0118 y0118 z0118 M  
 C-CT--0.066 0 x0119 y0119 z0119 M  
 H-H-0.109 0 x0120 y0120 z0120 M  
 H-H1-0.087 0 x0121 y0121 z0121 M  
 H-HC-0.019 0 x0122 y0122 z0122 M  
 H-HC-0.024 0 x0123 y0123 z0123 M  
 H-HC-0.024 0 x0124 y0124 z0124 M  
 H-HC-0.088 0 x0125 y0125 z0125 M  
 H-HC-0.088 0 x0126 y0126 z0126 M  
 H-HC-0.088 0 x0127 y0127 z0127 M  
 H-HC-0.019 0 x0128 y0128 z0128 M  
 H-HC-0.019 0 x0129 y0129 z0129 M  
 H-HC-0.019 0 x0130 y0130 z0130 M  
 N-N--0.376 0 x0131 y0131 z0131 M  
 C-CT--0.057 0 x0132 y0132 z0132 M  
 C-C-0.383 0 x0133 y0133 z0133 M  
 O-O--0.450 0 x0134 y0134 z0134 M  
 C-CT-0.298 0 x0135 y0135 z0135 M  
 C-CT--0.319 0 x0136 y0136 z0136 M  
 C-CT--0.319 0 x0137 y0137 z0137 M  
 H-H-0.272 0 x0138 y0138 z0138 M  
 H-H1-0.042 0 x0139 y0139 z0139 M  
 H-HC--0.030 0 x0140 y0140 z0140 M  
 H-HC-0.079 0 x0141 y0141 z0141 M  
 H-HC-0.079 0 x0142 y0142 z0142 M  
 H-HC-0.079 0 x0143 y0143 z0143 M  
 H-HC-0.079 0 x0144 y0144 z0144 M  
 H-HC-0.079 0 x0145 y0145 z0145 M  
 H-HC-0.079 0 x0146 y0146 z0146 M  
 N-N--0.318 0 x0147 y0147 z0147 M  
 C-CT--0.087 0 x0148 y0148 z0148 M  
 C-C-0.597 0 x0149 y0149 z0149 M  
 O-O--0.568 0 x0150 y0150 z0150 M  
 C-CT-0.298 0 x0151 y0151 z0151 M  
 C-CT--0.319 0 x0152 y0152 z0152 M  
 C-CT--0.319 0 x0153 y0153 z0153 M  
 H-H-0.196 0 x0154 y0154 z0154 M  
 H-H1-0.085 0 x0155 y0155 z0155 M  
 H-HC--0.030 0 x0156 y0156 z0156 M  
 H-HC-0.079 0 x0157 y0157 z0157 M  
 H-HC-0.079 0 x0158 y0158 z0158 M  
 H-HC-0.079 0 x0159 y0159 z0159 M  
 H-HC-0.079 0 x0160 y0160 z0160 M  
 H-HC-0.079 0 x0161 y0161 z0161 M  
 H-HC-0.079 0 x0162 y0162 z0162 M  
 N-N--0.416 0 x0163 y0163 z0163 M  
 C-CT--0.060 0 x0164 y0164 z0164 M  
 C-C-0.421 0 x0165 y0165 z0165 M  
 O-O--0.568 0 x0166 y0166 z0166 M  
 C-CT-0.130 0 x0167 y0167 z0167 M  
 C-CT--0.043 0 x0168 y0168 z0168 M  
 C-CT--0.320 0 x0169 y0169 z0169 M  
 C-CT--0.066 0 x0170 y0170 z0170 M  
 H-H-0.272 0 x0171 y0171 z0171 M  
 H-H1-0.057 0 x0172 y0172 z0172 M  
 H-HC-0.019 0 x0173 y0173 z0173 M  
 H-HC-0.024 0 x0174 y0174 z0174 M  
 H-HC-0.024 0 x0175 y0175 z0175 M  
 H-HC-0.056 0 x0176 y0176 z0176 M  
 H-HC-0.088 0 x0177 y0177 z0177 M  
 H-HC-0.088 0 x0178 y0178 z0178 M  
 H-HC-0.019 0 x0179 y0179 z0179 M  
 H-HC-0.019 0 x0180 y0180 z0180 M  
 H-HC-0.019 0 x0181 y0181 z0181 M  
 N-N--0.416 0 x0182 y0182 z0182 M  
 C-CT--0.025 0 x0183 y0183 z0183 M  
 C-C-0.597 0 x0184 y0184 z0184 M  
 O-O--0.568 0 x0185 y0185 z0185 M  
 H-H-0.272 0 x0186 y0186 z0186 M  
 H-H1-0.070 0 x0187 y0187 z0187 M  
 H-H1-0.070 0 x0188 y0188 z0188 M  
 N-N--0.416 0 x0189 y0189 z0189 M  
 C-CT--0.025 0 x0190 y0190 z0190 M  
 C-C-0.597 0 x0191 y0191 z0191 M  
 O-O--0.568 0 x0192 y0192 z0192 M  
 H-H-0.272 0 x0193 y0193 z0193 M  
 H-H1-0.070 0 x0194 y0194 z0194 M  
 H-H1-0.070 0 x0195 y0195 z0195 M  
 N-N--0.416 0 x0196 y0196 z0196 M  
 C-CT--0.087 0 x0197 y0197 z0197 M  
 C-C-0.597 0 x0198 y0198 z0198 M  
 O-O--0.568 0 x0199 y0199 z0199 M  
 C-CT-0.219 0 x0200 y0200 z0200 M  
 C-CT--0.134 0 x0201 y0201 z0201 M  
 C-CT--0.088 0 x0202 y0202 z0202 M  
 H-H-0.272 0 x0203 y0203 z0203 M  
 H-H1-0.097 0 x0204 y0204 z0204 M  
 H-HC--0.018 0 x0205 y0205 z0205 M  
 H-HC-0.014 0 x0206 y0206 z0206 M  
 H-HC-0.025 0 x0207 y0207 z0207 M  
 H-HC-0.057 0 x0208 y0208 z0208 M  
 H-HC-0.011 0 x0209 y0209 z0209 M  
 H-HC-0.028 0 x0210 y0210 z0210 M  
 H-HC-0.016 0 x0211 y0211 z0211 M  
 C-CT--0.004 0 x0212 y0212 z0212 M  
 O-O--0.142 0 x0213 y0213 z0213 M  
 H-N3-0.025 0 x0214 y0214 z0214 M  
 C-C-0.245 0 x0215 y0215 z0215 M  
 H-CT-0.004 0 x0216 y0216 z0216 M  
 H-H1-0.022 0 x0217 y0217 z0217 M  
 C-CT--0.025 0 x0218 y0218 z0218 M  
 H-H-0.272 0 x0219 y0219 z0219 M  
 N-N--0.416 0 x0220 y0220 z0220 M  
 H-C-0.597 0 x0221 y0221 z0221 M  
 H-H1-0.045 0 x0222 y0222 z0222 M  
 H-H1-0.070 0 x0223 y0223 z0223 M  
 C-CT--0.087 0 x0224 y0224 z0224 M  
 O-O--0.568 0 x0225 y0225 z0225 M  
 H-N--0.416 0 x0226 y0226 z0226 M  
 C-C-0.597 0 x0227 y0227 z0227 M  
 H-CT-0.298 0 x0228 y0228 z0228 M  
 H-H1-0.097 0 x0229 y0229 z0229 M  
 C-CT-0.034 0 x0230 y0230 z0230 M  
 H-H-0.272 0 x0231 y0231 z0231 M  
 N-N--0.416 0 x0232 y0232 z0232 M  
 H-C-0.412 0 x0233 y0233 z0233 M  
 H-CT--0.182 0 x0234 y0234 z0234 M  
 H-H1-0.082 0 x0235 y0235 z0235 M  
 C-CT--0.024 0 x0236 y0236 z0236 M  
 O-O--0.353 0 x0237 y0237 z0237 M  
 H-N3-0.121 0 x0238 y0238 z0238 M  
 C-C-0.597 0 x0239 y0239 z0239 M  
 H-CT-0.031 0 x0240 y0240 z0240 M  
 H-H1-0.088 0 x0241 y0241 z0241 M  
 C-CT--0.087 0 x0242 y0242 z0242    
 H-H-0.272 0 x0243 y0243 z0243    
 N-N--0.416 0 x0244 y0244 z0244    
 H-C-0.597 0 x0245 y0245 z0245    
 H-CT-0.298 0 x0246 y0246 z0246    
 H-H1-0.097 0 x0247 y0247 z0247    
 C-CT--0.017 0 x0248 y0248 z0248 M  
 O-O--0.125 0 x0249 y0249 z0249 M  
 H-N--0.088 0 x0250 y0250 z0250 M  
 C-C-0.169 0 x0251 y0251 z0251 M  
 H-CT-0.028 0 x0252 y0252 z0252 M  
 H-H1-0.027 0 x0253 y0253 z0253 M  
 C-CT-0.021 0 x0254 y0254 z0254 M  
 H-H-0.272 0 x0255 y0255 z0255 M  
 N-N--0.416 0 x0256 y0256 z0256 M  
 H-C-0.235 0 x0257 y0257 z0257 M  
 H-CT--0.045 0 x0258 y0258 z0258 M  
 H-H1-0.082 0 x0259 y0259 z0259 M  
 C-CT--0.019 0 x0260 y0260 z0260 M  
 O-O--0.212 0 x0261 y0261 z0261 M  
 H-N3-0.175 0 x0262 y0262 z0262 M  
 C-C-0.366 0 x0263 y0263 z0263 M  
 H-CT-0.012 0 x0264 y0264 z0264 M  
 H-H1-0.088 0 x0265 y0265 z0265 M  
 C-CT--0.007 0 x0266 y0266 z0266 M  
 H-H-0.144 0 x0267 y0267 z0267 M  
 N-N--0.188 0 x0268 y0268 z0268 M  
 H-C-0.295 0 x0269 y0269 z0269 M  
 H-H1-0.009 0 x0270 y0270 z0270 M  
 H-H1-0.032 0 x0271 y0271 z0271 M  
 C-CT--0.005 0 x0272 y0272 z0272 M  
 O-O--0.112 0 x0273 y0273 z0273 M  
 H-N--0.099 0 x0274 y0274 z0274 M  
 C-C-0.164 0 x0275 y0275 z0275 M  
 H-H1-0.006 0 x0276 y0276 z0276 M  
 H-H1-0.018 0 x0277 y0277 z0277 M  
 C-CT-0.003 0 x0278 y0278 z0278 M  
 H-H-0.044 0 x0279 y0279 z0279 M  
 N-N--0.086 0 x0280 y0280 z0280 M  
 H-C-0.028 0 x0281 y0281 z0281 M  
 H-CT--0.015 0 x0282 y0282 z0282 M  
 H-H1-0.012 0 x0283 y0283 z0283 M  
 C-CT--0.087 0 x0284 y0284 z0284 M  
 O-O--0.568 0 x0285 y0285 z0285 M  
 H-N--0.416 0 x0286 y0286 z0286 M  
 C-C-0.597 0 x0287 y0287 z0287 M  
 H-CT-0.298 0 x0288 y0288 z0288 M  
 H-H1-0.097 0 x0289 y0289 z0289 M  
 C-CT--0.087 0 x0290 y0290 z0290 M  
 H-H-0.272 0 x0291 y0291 z0291 M  
 N-N--0.416 0 x0292 y0292 z0292 M  
 H-C-0.597 0 x0293 y0293 z0293 M  
 H-CT-0.298 0 x0294 y0294 z0294 M  
 H-H1-0.097 0 x0295 y0295 z0295 M  
 Variables:
x0008 =     0.9499
y0008 =     3.3181
z0008 =     5.8014
x0009 =     3.5400
y0009 =     2.8072
z0009 =     7.2624
x0010 =     2.6720
y0010 =     3.6338
z0010 =     9.3784
x0011 =     1.9114
y0011 =     6.1403
z0011 =     8.1287
x0012 =     3.0789
y0012 =     5.3489
z0012 =     7.0159
x0013 =     3.5425
y0013 =     5.6747
z0013 =     8.7208
x0014 =     0.3394
y0014 =     5.1737
z0014 =     8.6812
x0015 =     0.4526
y0015 =     3.7825
z0015 =     9.8125
x0016 =     0.0710
y0016 =     3.5023
z0016 =     8.0793
 Constants:
x0001 =     1.9150
y0001 =     3.4090
z0001 =     6.0850
x0002 =     2.4510
y0002 =     2.9370
z0002 =     7.3480
x0003 =     1.8360
y0003 =     1.6080
z0003 =     7.6930
x0004 =     0.6770
y0004 =     1.3780
z0004 =     7.3840
x0005 =     2.1630
y0005 =     3.9910
z0005 =     8.4710
x0006 =     2.7150
y0006 =     5.3930
z0006 =     8.0530
x0007 =     0.6450
y0007 =     4.1220
z0007 =     8.7840
x0017 =     0.1850
y0017 =    10.8110
z0017 =     9.0070
x0018 =    -0.3470
y0018 =     9.4920
z0018 =     8.9950
x0019 =     0.3510
y0019 =     8.9510
z0019 =    10.2200
x0020 =     1.5840
y0020 =     8.9700
z0020 =    10.2830
x0021 =     0.1270
y0021 =     8.6880
z0021 =     7.7230
x0022 =    -0.4860
y0022 =     7.2310
z0022 =     7.6640
x0023 =    -0.1320
y0023 =     9.4110
z0023 =     6.4320
x0024 =     1.0985
y0024 =    11.0255
z0024 =     9.3807
x0025 =    -1.4454
y0025 =     9.4345
z0025 =     8.9821
x0026 =     1.2174
y0026 =     8.5995
z0026 =     7.8377
x0027 =    -0.6429
y0027 =     6.8581
z0027 =     8.6869
x0028 =     0.2076
y0028 =     6.5632
z0028 =     7.1321
x0029 =    -1.4484
y0029 =     7.2599
z0029 =     7.1321
x0030 =     0.8187
y0030 =     9.7860
z0030 =     6.0252
x0031 =    -0.8116
y0031 =    10.2562
z0031 =     6.6157
x0032 =    -0.5925
y0032 =     8.7195
z0032 =     5.7111
x0033 =     2.1980
y0033 =     8.3200
z0033 =     0.7710
x0034 =     2.7830
y0034 =     7.8530
z0034 =     2.0440
x0035 =     2.0050
y0035 =     6.7250
z0035 =     2.7130
x0036 =     0.7990
y0036 =     6.6320
z0036 =     2.5430
x0037 =     1.2362
y0037 =     8.1738
z0037 =     0.4996
x0038 =     2.8213
y0038 =     8.7048
z0038 =     2.7389
x0039 =     3.8041
y0039 =     7.4962
z0039 =     1.8439
x0040 =     2.6790
y0040 =     5.8770
z0040 =     3.4880
x0041 =     2.0050
y0041 =     4.6830
z0041 =     4.0380
x0042 =     2.7060
y0042 =     4.0060
z0042 =     5.2120
x0043 =     3.9310
y0043 =     4.0070
z0043 =     5.3390
x0044 =     3.6447
y0044 =     6.1249
z0044 =     3.6492
x0045 =     1.9136
y0045 =     3.9442
z0045 =     3.2282
x0046 =     1.0022
y0046 =     4.9862
z0046 =     4.3733
x0047 =     2.6170
y0047 =     0.7430
z0047 =     8.3240
x0048 =     2.1700
y0048 =    -0.5480
z0048 =     8.8500
x0049 =     2.7590
y0049 =    -0.8600
z0049 =    10.2240
x0050 =     3.9570
y0050 =    -0.6360
z0050 =    10.4760
x0051 =     2.5720
y0051 =    -1.6990
z0051 =     7.8970
x0052 =     4.0590
y0052 =    -1.8200
z0052 =     7.8390
x0053 =     2.0190
y0053 =    -3.0260
z0053 =     8.3340
x0054 =     3.5681
y0054 =     1.0729
z0054 =     8.4053
x0055 =     1.0761
y0055 =    -0.4718
z0055 =     8.9368
x0056 =     2.1537
y0056 =    -1.4486
z0056 =     6.9109
x0057 =     4.5139
y0057 =    -1.0598
z0057 =     8.4910
x0058 =     4.3987
y0058 =    -1.6669
z0058 =     6.8040
x0059 =     4.3583
y0059 =    -2.8223
z0059 =     8.1793
x0060 =     0.9212
y0060 =    -3.0048
z0060 =     8.2682
x0061 =     2.3208
y0061 =    -3.2242
z0061 =     9.3731
x0062 =     2.4108
y0062 =    -3.8192
z0062 =     7.6802
x0063 =     1.8990
y0063 =    -1.3840
z0063 =    11.1140
x0064 =     2.2910
y0064 =    -2.0300
z0064 =    12.4000
x0065 =     1.7300
y0065 =    -3.4450
z0065 =    12.4590
x0066 =     0.5170
y0066 =    -3.6570
z0066 =    12.4520
x0067 =     1.7560
y0067 =    -1.3000
z0067 =    13.6290
x0068 =     2.1440
y0068 =     0.1860
z0068 =    13.6130
x0069 =     2.2700
y0069 =    -2.0380
z0069 =    14.8890
x0070 =     1.5720
y0070 =     1.0030
z0070 =    14.7320
x0071 =     0.9402
y0071 =    -1.2906
z0071 =    10.8107
x0072 =     3.3906
y0072 =    -2.0104
z0072 =    12.4204
x0073 =     0.6561
y0073 =    -1.3147
z0073 =    13.6310
x0074 =     1.7971
y0074 =     0.6187
z0074 =    12.6630
x0075 =     3.2406
y0075 =     0.2527
z0075 =    13.6682
x0076 =     2.3993
y0076 =    -1.3162
z0076 =    15.7089
x0077 =     3.2350
y0077 =    -2.5161
z0077 =    14.6651
x0078 =     1.5404
y0078 =    -2.8054
z0078 =    15.1870
x0079 =     1.2263
y0079 =     0.3342
z0079 =    15.5340
x0080 =     0.7243
y0080 =     1.5952
z0080 =    14.3568
x0081 =     2.3462
y0081 =     1.6783
z0081 =    15.1252
x0082 =    -2.1560
y0082 =    -2.5460
z0082 =    12.1250
x0083 =    -2.6900
y0083 =    -1.2810
z0083 =    11.7490
x0084 =    -1.9340
y0084 =    -0.8440
z0084 =    10.4470
x0085 =    -0.7010
y0085 =    -1.0090
z0085 =    10.3210
x0086 =    -2.4690
y0086 =    -0.2920
z0086 =    12.9080
x0087 =    -2.9720
y0087 =     1.1210
z0087 =    12.5760
x0088 =    -3.1000
y0088 =    -0.8160
z0088 =    14.1700
x0089 =    -1.1743
y0089 =    -2.7065
z0089 =    12.3001
x0090 =    -3.7709
y0090 =    -1.3162
z0090 =    11.5477
x0091 =    -1.3840
y0091 =    -0.2058
z0091 =    13.0670
x0092 =    -3.7776
y0092 =     1.3960
z0092 =    13.2727
x0093 =    -3.3555
y0093 =     1.1404
z0093 =    11.5452
x0094 =    -2.1426
y0094 =     1.8370
z0094 =    12.6727
x0095 =    -2.8514
y0095 =    -0.1484
z0095 =    15.0082
x0096 =    -2.7171
y0096 =    -1.8258
z0096 =    14.3788
x0097 =    -4.1920
y0097 =    -0.8562
z0097 =    14.0440
x0098 =    -2.7020
y0098 =    -0.3380
z0098 =     9.4720
x0099 =    -2.1860
y0099 =     0.2550
z0099 =     8.2550
x0100 =    -2.8650
y0100 =     1.5690
z0100 =     7.9020
x0101 =    -3.9910
y0101 =     1.8650
z0101 =     8.3320
x0102 =    -3.6923
y0102 =    -0.4051
z0102 =     9.6587
x0103 =    -2.3401
y0103 =    -0.4529
z0103 =     7.4272
x0104 =    -1.1094
y0104 =     0.4402
z0104 =     8.3842
x0105 =    -2.1880
y0105 =     2.3830
z0105 =     7.1170
x0106 =    -2.8720
y0106 =     3.5230
z0106 =     6.5160
x0107 =    -2.1180
y0107 =     3.9900
z0107 =     5.3110
x0108 =    -0.8870
y0108 =     3.7900
z0108 =     5.2250
x0109 =    -1.2142
y0109 =     2.1497
z0109 =     6.9853
x0110 =    -2.9311
y0110 =     4.3410
z0110 =     7.2491
x0111 =    -3.8874
y0111 =     3.2233
z0111 =     6.2175
x0112 =     4.2630
y0112 =    -1.4240
z0112 =     0.8350
x0113 =     4.6760
y0113 =    -1.9670
z0113 =     2.1570
x0114 =     4.1650
y0114 =    -3.3900
z0114 =     2.2420
x0115 =     3.0210
y0115 =    -3.6610
z0115 =     1.9120
x0116 =     4.0590
y0116 =    -1.2060
z0116 =     3.3400
x0117 =     4.3790
y0117 =     0.2810
z0117 =     3.2860
x0118 =     4.5680
y0118 =    -1.8050
z0118 =     4.6740
x0119 =     3.9390
y0119 =     1.0210
z0119 =     4.5040
x0120 =     3.3182
y0120 =    -1.4625
z0120 =     0.4800
x0121 =     5.7706
y0121 =    -1.8817
z0121 =     2.2241
x0122 =     2.9664
y0122 =    -1.3160
z0122 =     3.2751
x0123 =     3.8734
y0123 =     0.7158
z0123 =     2.4112
x0124 =     5.4672
y0124 =     0.4007
z0124 =     3.1783
x0125 =     5.3663
y0125 =    -1.1669
z0125 =     5.0809
x0126 =     4.9618
y0126 =    -2.8163
z0126 =     4.4946
x0127 =     3.7370
y0127 =    -1.8567
z0127 =     5.3928
x0128 =     4.4121
y0128 =     0.5777
z0128 =     5.3926
x0129 =     2.8451
y0129 =     0.9535
z0129 =     4.5982
x0130 =     4.2357
y0130 =     2.0767
z0130 =     4.4181
x0131 =     0.3030
y0131 =    -2.6960
z0131 =     2.0750
x0132 =    -0.2530
y0132 =    -1.4370
z0132 =     1.7140
x0133 =     0.4990
y0133 =    -0.8860
z0133 =     0.4540
x0134 =     1.7230
y0134 =    -0.8090
z0134 =     0.4420
x0135 =    -0.0860
y0135 =    -0.4400
z0135 =     2.9440
x0136 =    -0.5620
y0136 =     1.0120
z0136 =     2.6400
x0137 =    -0.6890
y0137 =    -0.9980
z0137 =     4.1960
x0138 =     1.2898
y0138 =    -2.8727
z0138 =     2.1983
x0139 =    -1.3209
y0139 =    -1.5387
z0139 =     1.4705
x0140 =     0.9963
y0140 =    -0.3510
z0140 =     3.1192
x0141 =    -1.3557
y0141 =     1.2929
z0141 =     3.3479
x0142 =    -0.9517
y0142 =     1.0620
z0142 =     1.6126
x0143 =     0.2859
y0143 =     1.7049
z0143 =     2.7452
x0144 =     0.0433
y0144 =    -1.6484
z0144 =     4.6968
x0145 =    -1.5856
y0145 =    -1.5823
z0145 =     3.9417
x0146 =    -0.9670
y0146 =    -0.1726
z0146 =     4.8679
x0147 =     4.2780
y0147 =     4.4850
z0147 =    14.6280
x0148 =     4.8620
y0148 =     4.9660
z0148 =    13.3390
x0149 =     4.2710
y0149 =     6.2880
z0149 =    12.8850
x0150 =     3.1040
y0150 =     6.5160
z0150 =    13.0800
x0151 =     4.6930
y0151 =     3.8540
z0151 =    12.2180
x0152 =     3.2800
y0152 =     3.4020
z0152 =    12.0980
x0153 =     5.2220
y0153 =     4.2970
z0153 =    10.8270
x0154 =     3.3044
y0154 =     4.5575
z0154 =    14.8868
x0155 =     5.9330
y0155 =     5.1484
z0155 =    13.5113
x0156 =     5.3144
y0156 =     3.0103
z0156 =    12.5527
x0157 =     2.7991
y0157 =     3.4373
z0157 =    13.0867
x0158 =     3.2566
y0158 =     2.3711
z0158 =    11.7149
x0159 =     2.7415
y0159 =     4.0640
z0159 =    11.4039
x0160 =     6.1640
y0160 =     4.8509
z0160 =    10.9529
x0161 =     4.4773
y0161 =     4.9447
z0161 =    10.3413
x0162 =     5.3990
y0162 =     3.4088
z0162 =    10.2027
x0163 =     4.2680
y0163 =     9.5750
z0163 =     9.5490
x0164 =     4.6370
y0164 =    10.5900
z0164 =     8.5540
x0165 =     3.4750
y0165 =    11.5230
z0165 =     8.8230
x0166 =     2.3530
y0166 =    11.0270
z0166 =     8.7910
x0167 =     4.4870
y0167 =    10.0340
z0167 =     7.1170
x0168 =     4.8860
y0168 =     8.5530
z0168 =     7.0580
x0169 =     5.2560
y0169 =    10.9200
z0169 =     6.0910
x0170 =     4.4160
y0170 =     7.7620
z0170 =     5.8560
x0171 =     3.3150
y0171 =     9.2717
z0171 =     9.6902
x0172 =     5.6602
y0172 =    10.9883
z0172 =     8.6210
x0173 =     3.4266
y0173 =    10.0804
z0173 =     6.8283
x0174 =     4.4804
y0174 =     8.0641
z0174 =     7.9560
x0175 =     5.9847
y0175 =     8.5044
z0175 =     7.0781
x0176 =     6.0135
y0176 =    11.5179
z0176 =     6.6189
x0177 =     4.5475
y0177 =    11.5907
z0177 =     5.5829
x0178 =     5.7490
y0178 =    10.2758
z0178 =     5.3481
x0179 =     4.9964
y0179 =     6.8306
z0179 =     5.7809
x0180 =     4.5618
y0180 =     8.3600
z0180 =     4.9443
x0181 =     3.3488
y0181 =     7.5214
z0181 =     5.9705
x0182 =     6.7680
y0182 =    -0.3380
z0182 =     9.4720
x0183 =     7.2840
y0183 =     0.2550
z0183 =     8.2550
x0184 =     6.6050
y0184 =     1.5690
z0184 =     7.9020
x0185 =     5.4790
y0185 =     1.8650
z0185 =     8.3320
x0186 =     5.7777
y0186 =    -0.4051
z0186 =     9.6587
x0187 =     7.1299
y0187 =    -0.4529
z0187 =     7.4272
x0188 =     8.3606
y0188 =     0.4402
z0188 =     8.3842
x0189 =     7.2820
y0189 =     2.3830
z0189 =     7.1170
x0190 =     6.5980
y0190 =     3.5230
z0190 =     6.5160
x0191 =     7.3520
y0191 =     3.9900
z0191 =     5.3110
x0192 =     8.5830
y0192 =     3.7900
z0192 =     5.2250
x0193 =     8.2558
y0193 =     2.1497
z0193 =     6.9853
x0194 =     6.5389
y0194 =     4.3410
z0194 =     7.2491
x0195 =     5.5826
y0195 =     3.2233
z0195 =     6.2175
x0196 =     5.0560
y0196 =     7.1610
z0196 =    12.2700
x0197 =     4.4810
y0197 =     8.3710
z0197 =    11.6340
x0198 =     5.1340
y0198 =     8.9530
z0198 =    10.3610
x0199 =     6.3600
y0199 =     8.8700
z0199 =    10.1430
x0200 =     4.3210
y0200 =     9.6130
z0200 =    12.6310
x0201 =     3.4720
y0201 =     9.2470
z0201 =    13.7550
x0202 =     5.6630
y0202 =    10.1430
z0202 =    13.1310
x0203 =     6.0362
y0203 =     6.9181
z0203 =    12.2851
x0204 =     3.5262
y0204 =     7.9184
z0204 =    11.3282
x0205 =     3.8483
y0205 =    10.4263
z0205 =    12.0609
x0206 =     3.1241
y0206 =    10.1585
z0206 =    14.2631
x0207 =     4.0458
y0207 =     8.6291
z0207 =    14.4613
x0208 =     2.6052
y0208 =     8.6772
z0208 =    13.3889
x0209 =     5.5547
y0209 =    11.1985
z0209 =    13.4211
x0210 =     6.4119
y0210 =    10.0573
z0210 =    12.3299
x0211 =     5.9877
y0211 =     9.5544
z0211 =    14.0017
x0212 =     0.3390
y0212 =    13.0740
z0212 =     8.4090
x0213 =    -1.6500
y0213 =    11.8510
z0213 =     8.2300
x0214 =    -0.0280
y0214 =    13.6267
z0214 =     7.5316
x0215 =    -0.4710
y0215 =    11.8510
z0215 =     8.5330
x0216 =     0.2410
y0216 =    13.6333
z0216 =     9.3511
x0217 =     1.4095
y0217 =    12.8757
z0217 =     8.2516
x0218 =     0.1590
y0218 =     7.9090
z0218 =    12.4230
x0219 =    -1.4146
y0219 =     8.5333
z0219 =    11.0055
x0220 =    -0.4240
y0220 =     8.4720
z0220 =    11.1930
x0221 =    -0.3168
y0221 =     6.9769
z0221 =    12.7619
x0222 =     0.0622
y0222 =     8.6566
z0222 =    13.2240
x0223 =     1.2239
y0223 =     7.7050
z0223 =    12.2378
x0224 =     2.2350
y0224 =     9.4710
z0224 =    -1.4030
x0225 =     4.1500
y0225 =     9.1980
z0225 =     0.0210
x0226 =     2.5963
y0226 =    10.5070
z0226 =    -1.4815
x0227 =     2.9350
y0227 =     8.9880
z0227 =    -0.1290
x0228 =     2.5410
y0228 =     8.8411
z0228 =    -2.2513
x0229 =     1.1365
y0229 =     9.4154
z0229 =    -1.3918
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
x0236 =    -2.1160
y0236 =    -4.8860
z0236 =    12.3220
x0237 =    -4.1290
y0237 =    -3.5680
z0237 =    12.3980
x0238 =    -2.2211
y0238 =    -5.3286
z0238 =    13.3235
x0239 =    -2.9000
y0239 =    -3.6020
z0239 =    12.2730
x0240 =    -2.4235
y0240 =    -5.5679
z0240 =    11.5155
x0241 =    -1.0485
y0241 =    -4.6843
z0241 =    12.1495
x0242 =    -2.1550
y0242 =     5.2040
z0242 =     3.2200
x0243 =    -3.8204
y0243 =     4.6750
z0243 =     4.5642
x0244 =    -2.8270
y0244 =     4.6290
z0244 =     4.3880
x0245 =    -2.6480
y0245 =     6.1139
z0245 =     2.8471
x0246 =    -2.1255
y0246 =     4.4702
z0246 =     2.4011
x0247 =    -1.1404
y0247 =     5.4720
z0247 =     3.5500
x0248 =     4.5500
y0248 =    -0.6090
z0248 =    -1.4000
x0249 =     6.2640
y0249 =    -0.4990
z0249 =     0.2790
x0250 =     4.9609
y0250 =     0.3199
z0250 =    -1.8222
x0251 =     5.1010
y0251 =    -0.8300
z0251 =    -0.0170
x0252 =     4.8325
y0252 =    -1.4979
z0252 =    -1.9832
x0253 =     3.4567
y0253 =    -0.4878
z0253 =    -1.4054
x0254 =     4.5060
y0254 =    -5.6570
z0254 =     2.6390
x0255 =     5.8753
y0255 =    -4.0038
z0255 =     3.0205
x0256 =     4.9740
y0256 =    -4.3080
z0256 =     2.6810
x0257 =     5.1720
y0257 =    -6.3564
z0257 =     3.1656
x0258 =     4.4015
y0258 =    -5.9849
z0258 =     1.5942
x0259 =     3.5329
y0259 =    -5.6691
z0259 =     3.1518
x0260 =     0.2040
y0260 =    -4.9640
z0260 =     2.8010
x0261 =    -1.6900
y0261 =    -3.6770
z0261 =     2.0570
x0262 =    -0.3459
y0262 =    -5.2255
z0262 =     3.7171
x0263 =    -0.4700
y0263 =    -3.7170
z0263 =     2.2700
x0264 =     0.1277
y0264 =    -5.8010
z0264 =     2.0914
x0265 =     1.2747
y0265 =    -4.7795
z0265 =     2.9732
x0266 =     0.3030
y0266 =     0.0730
z0266 =    -1.7930
x0267 =    -1.2394
y0267 =    -0.6648
z0267 =    -0.4541
x0268 =    -0.2500
y0268 =    -0.5090
z0268 =    -0.5840
x0269 =    -0.2129
y0269 =     0.9904
z0269 =    -2.1129
x0270 =     0.2287
y0270 =    -0.6685
z0270 =    -2.6021
x0271 =     1.3609
y0271 =     0.3136
z0271 =    -1.6114
x0272 =     4.3150
y0272 =     3.4470
z0272 =    16.8150
x0273 =     6.2740
y0273 =     3.7410
z0273 =    15.4470
x0274 =     4.7792
y0274 =     2.5517
z0274 =    17.2544
x0275 =     5.0370
y0275 =     3.9040
z0275 =    15.5600
x0276 =     4.3409
y0276 =     4.2600
z0276 =    17.5555
x0277 =     3.2703
y0277 =     3.2162
z0277 =    16.5595
x0278 =     4.8820
y0278 =    13.4940
z0278 =     9.4650
x0279 =     2.7496
y0279 =    13.2918
z0279 =     9.2306
x0280 =     3.6340
y0280 =    12.8070
z0280 =     9.1770
x0281 =     4.7044
y0281 =    14.1354
z0281 =    10.3408
x0282 =     5.1402
y0282 =    14.1639
z0282 =     8.6316
x0283 =     5.6842
y0283 =    12.7592
z0283 =     9.6280
x0284 =     6.7800
y0284 =    -1.2810
z0284 =    11.7490
x0285 =     8.7690
y0285 =    -1.0090
z0285 =    10.3210
x0286 =     7.1926
y0286 =    -2.2584
z0286 =    12.0395
x0287 =     7.5360
y0287 =    -0.8440
z0287 =    10.4470
x0288 =     6.9379
y0288 =    -0.5744
z0288 =    12.5771
x0289 =     5.6991
y0289 =    -1.3162
z0289 =    11.5477
x0290 =     7.3150
y0290 =     5.2040
z0290 =     3.2200
x0291 =     5.6496
y0291 =     4.6750
z0291 =     4.5642
x0292 =     6.6430
y0292 =     4.6290
z0292 =     4.3880
x0293 =     6.8220
y0293 =     6.1139
z0293 =     2.8471
x0294 =     7.3445
y0294 =     4.4702
z0294 =     2.4011
x0295 =     8.3296
y0295 =     5.4720
z0295 =     3.5500
 
HrmStr1   * * 331. 1.09
HrmBnd1   * * * 20. 90.
AmbTrs    * * * * 0 0 0 0 0.0 0.0 1.15 0.0 3.0
VDW       MG 1.09 0.25
VDW       NT 1.8240 0.17

   1-0295 0
6-31G*  
****
 
END
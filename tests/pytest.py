import matplotlib
matplotlib.use('AGG')
import matplotlib.pyplot as plt


fig=plt.figure(figsize=(12, 14))
ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=5)
MODIFIED_FRAMESHIFT = 0
MODIFIED_NON_FRAMESHIFT=0
NON_MODIFIED_NON_FRAMESHIFT=0
patches, texts, autotexts =ax1.pie([MODIFIED_FRAMESHIFT,\
                   MODIFIED_NON_FRAMESHIFT,\
                  NON_MODIFIED_NON_FRAMESHIFT],\
                    labels=['Frameshift mutation\n(%d reads)' %MODIFIED_FRAMESHIFT,\
                    'In-frame mutation\n(%d reads)' % MODIFIED_NON_FRAMESHIFT,\
                    'Noncoding mutation\n(%d reads)' %NON_MODIFIED_NON_FRAMESHIFT],\
                    explode=(0.0, 0.0, 0.0),\
                        colors=[(0.89019608,  0.29019608,  0.2, 0.8), (0.99215686,  0.73333333,  0.51764706, 0.8), (0.99607843,  0.90980392,  0.78431373, 0.8)],\
                         autopct='%1.2f%%')

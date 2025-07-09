import glob as glob
import os

files = glob.glob('../../diapl2/WorkingDir/*.fits')

with open('../../diapl2/WorkingDir/images.txt', 'w') as f:
    for line in files:
        f.write(f"{os.path.basename(line)}\n")

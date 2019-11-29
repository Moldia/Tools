# run imagej stitching for image sets in the input file, one set of images in one directory
# parallel processing possible, but not implemented (need to avoid overwriting of macro files)
# Xiaoyan, 2017


import subprocess
import os


imlistfile = r"K:\stitchingdir.txt"
Processes = []
NextProcess = 0
stitchingwork =[]
ijdir = 'C:/Users/Public/Fiji.app/'


def correctinput(string):
    # string = string.encode('unicode-escape').decode()   # un-escape escape characters
    string = string.replace('\\', '/')  # change every unescaped backslash to slash
    string = string.replace('//', '/')  # change original double backslash (one for escaping) to one slash
    return string


def readdirfile(file):
    with open(file, 'r') as f:
        alldirs = [str(correctinput(line.strip('\n'))) for line in f]
    return alldirs


def renameijoutput(imdir):
    filels = [i for i in os.listdir(imdir) if os.path.isfile(os.path.join(imdir, i))]
    for i in filels:
        os.rename(os.path.join(imdir, i), os.path.join(imdir, i+'.tif'))


def findtostitch(imdir):
    """ For a given parent directory, find images to stitch """
    childdir = next(os.walk(imdir))[1]
    if 'FS_tophat_stack' in childdir:
        try:
            os.mkdir(imdir + '/Stitched')
        except:
            pass
    configfiles = os.listdir(imdir + '/FS_tophat_stack')
    configfiles = [i for i in configfiles if ('TileConfiguration' in i and 'registered' not in i)]
    tostitch = []
    for i in configfiles:
        outdir = i.split('Configuration')[1]
        outdir = outdir.split('.txt')[0]
        outdir = outdir.translate(None, '_')
        try:
            os.makedirs(imdir + '/Stitched/' + outdir)
            outdir = '/' + outdir
            tostitch.append([imdir, i, outdir])
        except:
            pass
    return tostitch


def newij():
    """ Start a new ImageJ incident and run stitching macro """
    global Processes
    global NextProcess
    if NextProcess < len(stitchingwork):
        strtojoin = ("type=[Positions from file]",
                     "order=[Defined by TileConfiguration]",
                     "directory=" + stitchingwork[NextProcess][0] + '/FS_tophat_stack',
                     "layout_file=" + stitchingwork[NextProcess][1],
                     "fusion_method=[Linear Blending]",
                     "regression_threshold=0.30",
                     "max/avg_displacement_threshold=2.50",
                     "absolute_displacement_threshold=3.50",
                     "compute_overlap", "computation_parameters=[Save memory (but be slower)]",
                     "image_output=[Write to disk]",
                     "output_directory=" + stitchingwork[NextProcess][0] + '/Stitched' + stitchingwork[NextProcess][2])
        with open(ijdir + 'macros/macroStitching.ijm', 'w') as f:
            f.write('run("Grid/Collection stitching", "' + " ".join(strtojoin) + '");\n' +
                    'eval("script", "System.exit(0);");\n')
        # imagej exits too early before output if run in headless
        ijcmd = " ".join((ijdir + 'ImageJ-win64', '-ij2', '-Xincgc', '-macro', 'macroStitching.ijm'))
        ijprocess = subprocess.Popen(ijcmd)
        Processes.append(ijprocess)
        NextProcess += 1


def runningij():
    global Processes
    global NextProcess
    for p in range(len(Processes)-1, -1, -1):     # check the processes in reverse order
        if Processes[p].poll() is not None:     # if the process hasn't finished will return None
            del Processes[p]
            renameijoutput(stitchingwork[NextProcess-1][0] + '/Stitched' + stitchingwork[NextProcess-1][2])
    while (len(Processes) < 1) and (NextProcess < len(stitchingwork)):  # more to do and some spare slots
        print("%d out of %d" % (NextProcess+1, len(stitchingwork)))
        newij()


alldirs = readdirfile(correctinput(imlistfile))
for i in alldirs:
    tileconfigs = findtostitch(i)
    for j in (tileconfigs):
        stitchingwork.append(j)

print("Start stitching..")
runningij()
while len(Processes) > 0:
    runningij()
renameijoutput(stitchingwork[-1][0] + '/Stitched' + stitchingwork[-1][2])

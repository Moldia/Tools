# run imagej stitching for image sets in the input file, one set of images in one directory
# parallel processing possible
# Xiaoyan, 2017


import subprocess
import os


imlistfile = r"G:\DemoData\NatureMethods\Preprocessing\images_to_stitch.txt"
Processes = []
NextProcess = 0
stitchingwork =[]
ijdir = 'E:/Programs/Fiji.app/'


def correctinput(string):
    # string = string.encode('unicode-escape').decode()   # un-escape escape characters
    string = string.replace('\\', '/')  # change every unescaped backslash to slash
    string = string.replace('//', '/')  # change original double backslash (one for escaping) to one slash
    return string


def readdirfile(file):
    with open(file, 'r') as f:
        alldirs = [str(correctinput(line.strip('\n'))) for line in f]
    return alldirs


def renameijoutput(stitched):
    filels = [i for i in os.listdir(stitched[2]) if os.path.isfile(os.path.join(stitched[2], i)) and i[-4:] != '.tif']
    for i in filels:
        os.rename(os.path.join(stitched[2], i), os.path.join(stitched[2].split(stitched[1])[0], stitched[1] + "_c" + i[5] + '.tif'))


def findtostitch(imdir):
    """ For a given parent directory, find images to stitch """
    childdir = next(os.walk(imdir))[1]
    tostitch = []
    for i in childdir:
        if i != 'Stitched':
            try:
                os.makedirs(os.path.join(imdir, 'Stitched', i))
            except:
                pass
            tostitch.append([correctinput(os.path.join(imdir, i)), i, correctinput(os.path.join(imdir, 'Stitched', i))])
    return tostitch


def newij():
    """ Start a new ImageJ incident and run stitching macro """
    global Processes
    global NextProcess
    if NextProcess < len(stitchingwork):
        strtojoin = ("type=[Grid: snake by rows]",
                     "order=[Right & Down                ]",
                     "grid_size_x=4 grid_size_y=4",
                     "tile_overlap=10",
                     "first_file_index_i=1",
                     "directory=" + stitchingwork[NextProcess][0],
                     "file_names=tile{i}.tif",
                     "output_textfile_name=TileConfiguration.txt",
                     "fusion_method=[Max. Intensity]",
                     "regression_threshold=0.30",
                     "max/avg_displacement_threshold=2.50",
                     "absolute_displacement_threshold=3.50",
                     "compute_overlap", "computation_parameters=[Save memory (but be slower)]",
                     "image_output=[Write to disk]",
                     "output_directory=" + stitchingwork[NextProcess][2])
        with open(ijdir + 'macros/macroStitching' + str(NextProcess) + '.ijm', 'w') as f:
            f.write('run("Grid/Collection stitching", "' + " ".join(strtojoin) + '");\n' +
                    'eval("script", "System.exit(0);");\n')
        # imagej exits too early before output if run in headless
        ijcmd = " ".join((ijdir + 'ImageJ-win64', '--ij2', '-macro',
                          'macroStitching' + str(NextProcess) + '.ijm'))
        ijprocess = subprocess.Popen(ijcmd)
        Processes.append(ijprocess)
        print("%d out of %d" % (NextProcess+1, len(stitchingwork)))
        NextProcess += 1


def runningij():
    global Processes
    global NextProcess
    for p in range(len(Processes)-1, -1, -1):     # check the processes in reverse order
        if Processes[p].poll() is not None:     # if the process hasn't finished will return None
            del Processes[p]
    while (len(Processes) < 3) and (NextProcess < len(stitchingwork)):  # more to do and some spare slots
        newij()


alldirs = readdirfile(correctinput(imlistfile))
for i in alldirs:
    tostitch = findtostitch(i)
    for j in tostitch:
        stitchingwork.append(j)

print("Start stitching..")
newij()
while len(Processes) > 0:
    runningij()

for c, i in enumerate(stitchingwork):
    renameijoutput(i)
    os.remove(os.path.join(ijdir, 'macros', 'macroStitching' + str(c) + '.ijm'))

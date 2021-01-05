import os

def normalize_files(files):
    for file in files:
        if int(file.split("_")[file.count("_")].replace(".stl","")) < 10:
            os.rename(file, file[0:file.rindex("_")+1] + "00" + file[file.rindex(".")-1:])
        elif int(file.split("_")[file.count("_")].replace(".stl","")) < 100:
            os.rename(file, file[0:file.rindex("_")+1] + "0" + file[file.rindex(".")-2:])
        else:
            os.rename(file, file)

def rename_files(files, i):
    extension = ".stl"
    template = "updated_electrode"
    if files[len(files)//2][0:8] != "electrode":
        template = "electrode"

    for file in files:
        if i < 10:
            temp = "00" + str(i)
        elif i < 100:
            temp = "0" + str(i)
        else:
            temp = str(i)
        new_name = template + temp + extension
        os.rename(file, new_name)
        print("renamed " + file + " to " + new_name)
        i = i + 1

    return i

dir = "Sextupole"
script_path = os.path.realpath(__file__)
dir_path = os.path.dirname(script_path) + "\\" + dir
os.chdir(dir_path)

rename_files(os.listdir(),1)
# normalize_files(os.listdir())
# files = os.listdir()
# rename_files(files[0:8], rename_files(files[8:], 1))
# print(files)

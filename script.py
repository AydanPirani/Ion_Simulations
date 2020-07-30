import os

dir = "STL"

i = 1
extension = ".stl"
template = "updated_electrode"

script_path = os.path.realpath(__file__)
dir_path = os.path.dirname(script_path) + "\\" + dir
files = os.listdir(dir_path)

if files[0][0:9] != "electrode":
    template = "electrode"

os.chdir(dir_path)

for file in files:
    if file[-4:] == ".stl":
        new_name = template + str(i) + extension
        os.rename(file, new_name)
        print("renamed " + file)
    else:
        print("kept " + file)
    i = i + 1

print()
print("done! renamed " + str(i-1) + " files")

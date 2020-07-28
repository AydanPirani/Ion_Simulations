import os

dir = "STL"

i = 1
extension = ".stl"
template = "electrode"

script_path = os.path.realpath(__file__)
dir_path = os.path.dirname(script_path) + "\\" + dir
files = os.listdir(dir_path)

os.chdir(dir_path)

for file in files:
    new_name = template + str(i) + extension
    os.rename(file, new_name)
    print("renamed " + file)
    i = i + 1
print("Done!")

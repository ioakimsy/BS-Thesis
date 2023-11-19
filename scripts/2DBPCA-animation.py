import glob
import os
from PIL import Image as pilim

def convertGIFtoPNG(imageLoadRelativePath, gifSaveRelativePath):
    '''
    imageLoadRelativePath: str, relative path to folder where to load png images
    gifSaveRelativePath:   str, relative path to folder where to save gif
    
    Example:
        >>> imageLoadRelativePath = './insert/path/to/load/images'
        >>> gifSaveRelativePath   = './insert/path/to/save/'
        >>> convertGIFtoPNG(imageLoadRelativePath, gifSaveRelativePath)
    '''
    
    loadImagePATH = f'{imageLoadRelativePath}/*.png'
    saveGIFPATH   = f'{gifSaveRelativePath}/output-gif.gif'
    gif, *imgs    = [pilim.open(_file) for _file in sorted(glob.glob(loadImagePATH), key=os.path.getmtime)]
    gif.save(fp=saveGIFPATH, format='GIF', append_images=imgs, save_all=True, duration=100, loop=0)

# List of parameters
sizes = [128]
seat_configs = ["random"]
Λs = [0.25,0.5,0.75]
steady_state_tolerance = 10
n_trials = 3
n_learned = 4

for seat_config in seat_configs:
    for λ₀ in Λs:
        for class_size in sizes :
            for trial in 1:n_trials:
                println("$seat_config 	$λ₀ 	$class_size 	$trial")


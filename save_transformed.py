# Adds the command save_transformed
# Usage: save_transformed object, file
def save_transformed(object,file):
    m = cmd.get_view(0)
    ttt = [m[0], m[1], m[2], 0.0,
           m[3], m[4], m[5], 0.0,
           m[6], m[7], m[8], 0.0,
           0.0,   0.0,  0.0, 1.0]
    cmd.transform_object(object,ttt,transpose=1)
    cmd.save(file,object)

cmd.extend('save_transformed',save_transformed)
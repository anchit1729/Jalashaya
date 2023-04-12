//
// Created by Anchit Mishra on 2023-04-08.
//

#ifndef FLUIDSIM_3D_UTILS_H
#define FLUIDSIM_3D_UTILS_H

class Utils {
public:
    static float clamp(float val, float min, float max) {
        if (val < min) return min;
        if (val > max) return max;
        return val;
    }
};

#endif //FLUIDSIM_3D_UTILS_H

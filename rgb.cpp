//
// Created by dell on 4/2/2021.
//
#include "rgb.h"
RGB::RGB() {

}
RGB& RGB::operator=(RGB rgb) {
    this->r=rgb.r;
    this->g=rgb.g;
    this->b=rgb.b;
    return *this;
}

RGB::RGB(int r, int g, int b) {
    this->r=r;
    this->g=g;
    this->b=b;

}

bool RGB::operator==(RGB rgb) {
    if(this->r==rgb.r && this->g==rgb.g && this->b==rgb.b)
        return true;
    else
        return false;
}

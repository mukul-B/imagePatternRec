#ifndef RGB_H
#define RGB_H

// an example -- you would need to add more functions

struct RGB {
    RGB(int, int, int);
    RGB();
    RGB& operator=(RGB);
    bool operator==(RGB);
    int r, g, b;
} ;

#endif
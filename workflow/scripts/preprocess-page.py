# References
# https://tesseract-ocr.github.io/tessdoc/ImproveQuality.html
# https://nanonets.com/blog/ocr-with-tesseract/#preprocessingfortesseract

import cv2
import numpy as np


# get grayscale image
def get_grayscale(image):
    return cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)


# noise removal
def remove_noise(image):
    return cv2.medianBlur(image, 5)


# thresholding
def thresholding(image):
    return cv2.threshold(image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]


# opening - erosion followed by dilation
def opening(image):
    kernel = np.ones((5, 5), np.uint8)
    return cv2.morphologyEx(image, cv2.MORPH_OPEN, kernel)


# canny edge detection
def canny(image):
    return cv2.Canny(image, 100, 200)


# dilation
def dilate(image):
    kernel = np.ones((5, 5), np.uint8)
    return cv2.dilate(image, kernel, iterations=1)


# erosion
def erode(image):
    kernel = np.ones((5, 5), np.uint8)
    return cv2.erode(image, kernel, iterations=1)


# template matching
def match_template(image, template):
    return cv2.matchTemplate(image, template, cv2.TM_CCOEFF_NORMED)


if __name__ == "__main__":

    sys.stderr = open(snakemake.log[0], "w")

    image = cv2.imread(snakemake.input[0])

    # TODO Check which preprocessing techniques deliver the best results
    processed_image = get_grayscale(image)
    # processed_image = remove_noise(processed_image)
    processed_image = thresholding(processed_image)
    # processed_image = opening(processed_image)
    # processed_image = canny(processed_image)

    # TODO add deskewing
    # image = deskew(image)

    cv2.imwrite(snakemake.output[0], processed_image)

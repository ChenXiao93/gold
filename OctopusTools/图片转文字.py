import pytesseract
from PIL import Image

def perform_ocr(image_path):
    image = Image.open(image_path)

    # text = pytesseract.image_to_string(image, lang = 'eng')
    print(image)
    return

image_path = "./1.png"
perform_ocr(image_path)

# print(text)
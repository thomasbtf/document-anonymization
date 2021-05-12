import pytesseract
import sys
import pandas as pd
import io
import time
import Levenshtein
import scipy.cluster.hierarchy as hi_cl
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import cv2
import typing

start = time.perf_counter()

sens_data = {"name": "Karadag", "first_name": "Murat", "street": "Am alten Sportplatz", "city": "Marl", "birthDate": "07.07.1985"}


detected_text  = pytesseract.image_to_data(sys.argv[1], lang="deu")
detected_text_df = pd.read_csv(io.StringIO(detected_text), sep='\t', quoting=3, error_bad_lines=False) # quoting=csv.QUOTE_NONE
filled_text_df = detected_text_df.fillna("")
# print(filled_text_df)



detected_boxes  = pytesseract.image_to_boxes(sys.argv[1])
print("#################################################")
print("Detected Boxes")
print("#################################################")
print(detected_boxes)
detected_orient  = pytesseract.image_to_osd(sys.argv[1])
#print("#################################################")
#print("Detected Orientation") only good for detection of 0°,90°,180°,270° rotation of whole image
#print("#################################################")
#print(detected_orient)

results_container = []
results_df = filled_text_df.iloc[0:0,:].copy()


# def match_word(data_dict, row, target_df):
#     text = row["text"]
#     matches = {k:v for (k,v) in sens_data.items() if v in text}
#     if any(matches):
#         target_df = target_df.append(row)
#         results_container.append(target_df)


# filled_text_df.apply(lambda row: match_word(sens_data, row, results_df), axis=1)

def segment_dataframe(detected_text_df: pd.DataFrame, col_key: str) -> dict:
    # get number of blocks
    max_blocks = detected_text_df[col_key].max()
    #print(max_blocks)
    #print(type(max_blocks))
    #block_tainer = []

    # add individual blocks to dict
    block_dict = {}
    for b in range (max_blocks + 1):
        block_in = detected_text_df[col_key].eq(b)
        block_out = detected_text_df[block_in]
        block_filt = block_out.dropna()
        block_sort = block_filt.sort_values(by=["left", "top"])
        #print(type(block_filt))
        #block_tainer.append(block_filt)
        block_dict[b] = block_sort
        #print(block_sort)

    # drop blocks that don`t contain any matches with sensitive data
    def block_match(block_df, sens_data):
        match_count = []
        
        def row_match(sens_data, row):
            text = row["text"]
            matches = {k:v for (k,v) in sens_data.items() if v in text}
            if any(matches):
                match_count.append(1)
                #print("match")

        
        block_df.apply(lambda row: row_match(sens_data, row), axis=1)
        
        return match_count
        
    block_dict = {k:v for (k,v) in block_dict.items() if len(block_match(v, sens_data)) > 0}
    #print("BLOCK_DICT:")
    #print(block_dict)
    return block_dict

block_dict = segment_dataframe(detected_text_df, "block_num")

def get_bounds(df: pd.DataFrame) -> tuple:
    # get block boundaries
#    block_bounds = {}
#    for block in block_dict:
#        top_min = block_dict[block].top.min()
#        top_max = block_dict[block].top.max()
#        left_min = block_dict[block].left.min()
#        left_max = block_dict[block].left.max()
#        boundaries = (top_min, top_max, left_min, left_max)
#        block_bounds[block] = boundaries
    top_min = df.top.min()
    top_max = df.top.max() + df.at[df.top.idxmax(), "height"]
    print("INDEX", df.top.idxmax())
    left_min = df.left.min()
    left_max = df.left.max() + df.at[df.left.idxmax(), "width"]
    print("INDEX", df.left.idxmax())
    boundaries = (top_min, top_max, left_min, left_max)
    # block_bounds[block] = boundaries
    return boundaries
    

# cluster data
for block in block_dict:
    left_vals = pd.DataFrame(block_dict[block].left)
    print(type(left_vals))
    print("left_vals Dataframe:")
    print(left_vals)
    # left_vals = np.array(left_vals)
    left_vals = left_vals.to_numpy()
    print("left_vals np_array:")
    print(left_vals)
    fig = plt.figure(figsize=(15, 15))
    ax1 = fig.add_axes([0.1, 0.1, 0.2, 0.6])
    #hi_cl.linkage(left_vals)
    link = hi_cl.linkage(left_vals, method='single')
    print(link)
    cut_tree = hi_cl.cut_tree(link, height=len(left_vals))
    print("cut_tree:")
    print(cut_tree)
    tree = hi_cl.dendrogram(link, orientation='left')
    fclust = hi_cl.fcluster(link, t=100, criterion='distance')
    print(fclust)
    print(type(fclust))
    subclust = pd.Series(fclust, index=block_dict[block].index)
    print(subclust)
    block_dict[block] = block_dict[block].assign(subclust=subclust)
    print("block:", block_dict[block])
    #plt.show()

clustainer = []
for block in block_dict:
    block_bounds = get_bounds(block_dict[block])
    print("block_bounds", block_bounds)
    sub_dict = segment_dataframe(block_dict[block], "subclust")
    clustainer.append(sub_dict)
    print(clustainer)
    
target_areas = []
for sub_dict in clustainer:
    for df in sub_dict:
        stamp_bounds = get_bounds(sub_dict[df])
        print("stamp_bounds", stamp_bounds)
        target_areas.append(stamp_bounds)
        


for key in sens_data:
    df_in = filled_text_df.text.str.contains(sens_data[key])
    df_out = filled_text_df[df_in]
    df = pd.DataFrame(df_out)
    results_container.append(df)


for df in results_container:
    results_df = results_df.append(df)
print(results_df)

# split blocks based on proximity clustering of text
# verify the tilt/angle/inclination of individual text snippets by relating image_to_text data to image_to_blocks data
# (maybe also verify font size through image_to_box data and relate to height of boxes from image_to_text data)
#calculate an average for text tilt over whole document
# group text with similar angles together
# using open_cv image segmentation to detect stamps/stickers

def redact(bounds: tuple, img: typing.Any) -> typing.Any:
   # (x, y, w, h) = (bounds[2], bounds[0], (bounds[3]-bounds[2]), bounds[1]-bounds[0])
   # (x, y, w, h) = (bounds[2], bounds[0], (bounds[3]-bounds[2]) + 100, bounds[1]-bounds[0] + 100)
    (x, y, w, h) = (bounds[2], bounds[0], (bounds[3]), bounds[1])
    # img = cv2.rectangle(img, (x, y), (x + w, y + h), (0, 0, 0), -1)
    img = cv2.rectangle(img, (x, y), (w, h), (0, 0, 0), -1)
    return img

    top_min = df.top.min()
    top_max = df.top.max()
    left_min = df.left.min()
    left_max = df.left.max()
    boundaries = (top_min, top_max, left_min, left_max)



img = cv2.imread(sys.argv[1])
for target in target_areas:
    img = redact(target, img)
cv2.imwrite(sys.argv[1] + "mod.tiff", img)
  
end = time.perf_counter()
print(end - start)
    

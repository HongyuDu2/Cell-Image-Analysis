import numpy as np
from PIL import Image
import os
import pandas as pd


def main():
    # directory of the project
    base_dir = r'C:\Users\duho\Desktop\FOVs\SLE39\colon\segmentation'
    # directory where deepcell output is stored
    seg_mask_dir = os.path.join(base_dir, 'deepcell_output')
    # directory where the cell table with 'cell_type' column is stored
    data_dir = r'C:\Users\duho\Desktop'
    # directory where the output phenotyping masks will be stored
    output_dir = os.path.join(base_dir, 'phenotyping_mask')

    # list all whole cell masks
    all_seg_masks = [mask_name for mask_name in os.listdir(seg_mask_dir) if mask_name.endswith('whole_cell.tiff')]


    # read the cell table with cell_type
    df = pd.read_csv(os.path.join(data_dir, 'cell_segmentation_table_new.csv'))
    # map every cell_type to an integer number and generate a new column 'cell_type_int'
    cell_type_to_numbers = {cell_type: i for i, cell_type in enumerate(df['cell_type'].unique())}
    df['cell_type_int'] = df['cell_type'].map(cell_type_to_numbers)

    # Write the dictionary to a .txt file
    with open(os.path.join(data_dir, 'cell_type_to_numbers.txt'), 'w') as f:
        for key, value in cell_type_to_numbers.items():
            f.write(f'{key}: {value}\n')



    # loop over every fov and generate a cellType mask
    for seg_mask in all_seg_masks:
        print(f'generating phenotyping mask for {seg_mask}')

        tiff_name = seg_mask.split('_whole_cell')[0]

        # generate the path to mask.tiff
        seg_mask_path = os.path.join(seg_mask_dir, seg_mask)

        # read the seg_mask
        seg_mask_array = read_image(seg_mask_path)

        # map cell_ids (labels) to cellTypes
        cellID_cellType_dict = generate_cell_id_cellType_map(df, tiff_name)

        # generate phenotype_mask
        pheno_mask_array = generate_pheno_mask(cellID_cellType_dict, seg_mask_array)

        # write the
        save_pheno_mask_to_image(output_dir, tiff_name, pheno_mask_array)

def read_image(image_path):
    image = Image.open(image_path)

    # Convert the image to a NumPy array for manipulation
    image_array = np.array(image)
    return image_array


def generate_cell_id_cellType_map(df, tiff_name):

    df = df[df.fov == tiff_name]
    required_cols = ['label', 'cell_type_int']
    # Check if all required columns are present in the DataFrame
    if all(col in df.columns for col in required_cols):
        df = df[required_cols]
        return dict(zip(df['label'], df['cell_type_int']))
    else:
        print('Required columns (label and cellType_int) are not in the df')


def generate_pheno_mask(cellID_cellType_dict, image_array):
    # Create a copy of the image_array to avoid modifying the original
    modified_image_array = np.copy(image_array)

    for i in range(modified_image_array.shape[0]):
        for j in range(modified_image_array.shape[1]):
            pixel_value = modified_image_array[i, j]

            if pixel_value != 0:
                try:
                    # Attempt to retrieve cell type from the dictionary
                    cell_type = cellID_cellType_dict[pixel_value]
                    modified_image_array[i, j] = cell_type
                except KeyError:
                    # Raise an error if pixel_value is not found in the dictionary
                    raise ValueError(f"Pixel value {pixel_value} not found in the dictionary")

    # Convert the NumPy array back to an image
    pheno_mask = Image.fromarray(modified_image_array)
    return pheno_mask
def save_pheno_mask_to_image(output_dir, tiff_name, pheno_mask_array):
    # Save the new image as a 32-bit grayscale TIFF
    new_image_path = os.path.join(output_dir, f'{tiff_name}_cell_type.tiff')
    pheno_mask_array.save(new_image_path, format='TIFF')
    print("New image saved successfully.")

if __name__ == '__main__':
    main()

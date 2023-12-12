def download_DMW(yyyymmddhhmn,Ch,path_dest):
    product_name = 'ABI-L2-DMWF'
    os.makedirs(path_dest, exist_ok=True)

    year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%Y')
    day_of_year = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%j')
    hour = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%H')
    min = datetime.strptime(yyyymmddhhmn, '%Y%m%d%H%M').strftime('%M')

    # AMAZON repository information
    # https://noaa-goes16.s3.amazonaws.com/index.html
    bucket_name = 'noaa-goes16'

    # Initializes the S3 client
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    # -----------------------------------------------------------------------------------------------------------
    # File structure
    prefix = f'{product_name}/{year}/{day_of_year}/{hour}/OR_{product_name}-M6C{int(Ch):02.0f}_G16_s{year}{day_of_year}{hour}{min}'

    # Seach for the file on the server
    s3_result = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    # -----------------------------------------------------------------------------------------------------------
    # Check if there are files available
    if 'Contents' not in s3_result:
        # There are no files
        print(f'No files found for the date: {yyyymmddhhmn}, Product-{product_name}')
        return -1
    else:
        # There are files
        for obj in s3_result['Contents']:
            key = obj['Key']
            # Print the file name
            file_name = key.split('/')[-1].split('.')[0]

            # Download the file
            if os.path.exists(f'{path_dest}/{file_name}.nc'):
                print(f'File {path_dest}/{file_name}.nc exists')
            else:
                print(f'Downloading file {path_dest}/{file_name}.nc')
                s3_client.download_file(bucket_name, key, f'{path_dest}/{file_name}.nc')
    return f'{file_name}.nc'

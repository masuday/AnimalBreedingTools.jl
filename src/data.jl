"""
    data = load_data_set(name)

Load a small data set for numerical testing.
"""
function load_data_set(dataset="")
   name = lowercase(dataset)
   if name=="small5ren" || name=="small5" || name=="small5ren2"
      if name=="small5ren"
         pedlist = [  10  16  12  16   8  10  15  16  18  15  18  14  10  18  14
                       7  17   1  17  11   7   9  17  19  11  19   4   5  19   2];
         genid = [ 14,  2,  9,  5,  7,  13,  1,  3,  6]
      elseif name=="small5ren2"
         pedlist = [  10  16  12  16   8  10  15  16  16  15  17  14  10  17  14
                       7  17   1  17  11   7   9  17  16  11  16   4   5  17   2];
         genid = [ 14,  2,  9,  5,  7,  13,  1,  3,  6]   
      else
         pedlist =[
            0  0  0  0  0  0  2  1  2  7  7  11  11   9  11
            0  0  0  0  0  0  5  4  3  6  4   8  10  13  10
         ]
         genid = [2,5,6,8,10,12,13,14,15]   
      end
      G = [
         1.19677     0.542311     0.150994    0.457376    0.162172    0.0959141  -0.0896063   0.513584    -0.0896064
         0.542311    1.45157      0.141889    0.475572    0.535284    0.0322102  -0.126011    0.00187444  -0.126011
         0.150994    0.141889     0.941952   -0.0249475   0.373563    0.0175262   0.410923    0.511502     0.410924
         0.457376    0.475572    -0.0249475   1.09059     0.22076     0.652166    0.0564771  -0.0585542    0.0564771
         0.162172    0.535284     0.373563    0.22076     0.778136   -0.0314916   0.318484    0.331571     0.318484
         0.0959141   0.0322102    0.0175262   0.652166   -0.0314916   1.34492     0.492344   -0.118302     0.492345
        -0.0896063  -0.126011     0.410923    0.0564771   0.318484    0.492344    0.798901    0.391289     0.77654
         0.513584    0.00187444   0.511502   -0.0585542   0.331571   -0.118302    0.391289    1.03575      0.380108
        -0.0896064  -0.126011     0.410924    0.0564771   0.318484    0.492345    0.77654     0.380108     0.798901
      ]
      y =[
         -1.671506
         0.467041
        -1.949929
         0.537337
        -0.888975
        -1.672478
         0.893532
        -1.278057
         0.040003
        -1.1635
         1.364332
        -0.475337
         0.596651
        -1.775682
         0.68465
      ]
      X = ones(length(y),1)
      A = AnimalBreedingDataSetRaw(pedlist,genid,G,y,X)
   else
      error("specify one of the following dataset:\n small5, small5ren, or small5ren2")
   end
   return A
end

# MeshModelingforCrimeDensity

Application of content-adaptive mesh modeling to crime spatial density estimation problems.

Reference: J. G. Brankov, Y. Yang, and M. N. Wernick, “Tomographic Image Reconstruction Based on a Content-Adaptive Mesh Model,” IEEE Trans. Med. Imaging, vol. 23, no. 2, pp. 202–212, 2004

1.Using mesh grid model to represent crime density

Assumed "true" density map
![alt text](https://github.com/xiaomuliu/MeshModelingforCrimeDensity/blob/master/TrueDen.png)

Contructed Mesh grid
![alt text](https://github.com/xiaomuliu/MeshModelingforCrimeDensity/blob/master/Mesh_1079.png)

Representation error of mesh grid with different nodes (MDL)
![alt text](https://github.com/xiaomuliu/MeshModelingforCrimeDensity/blob/master/Rep_MDL.png)

Using observed incidents to estimate density (comparing mesh, equivalent-sized pixel-grid, ML-estimation, MAP-estimation and different model configurations)
![alt text](https://github.com/xiaomuliu/MeshModelingforCrimeDensity/blob/master/Recon_RMSE.png)

Reconstructed density
![alt text](https://github.com/xiaomuliu/MeshModelingforCrimeDensity/blob/master/Recon_MAP10-3_node834.jpeg)

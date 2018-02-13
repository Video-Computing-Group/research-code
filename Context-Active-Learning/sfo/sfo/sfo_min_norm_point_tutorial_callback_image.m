% Helper callback function for visualizing inference in Ising model

function sfo_min_norm_point_tutorial_callback_image(Abest, img)
img(Abest)=1;
imshow(img)
pause(.1);

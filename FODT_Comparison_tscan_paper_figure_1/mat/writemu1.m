function data=writemu1(foo,var, data, dim1, dim2, dim3)
% function data=writemu1(foo,var, data, dim1, dim2, dim3)

foo1=[foo '.dat'];
fid=datOpen(foo1,'w+');
dims=[2 dim1 dim2 dim3];
write_float_array(fid, var, data, length(dims), dims);
fclose(fid);

return;


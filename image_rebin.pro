pro image_rebin, ima_in, err_in, dx, dy, ima_out, err_out
; This code aims to rebin an image with equal pixel interval 
; to another interval with the total flux unchanged. The dx and 
; dy could be different but themselves should be constant. 
; We assume the output image has the same start point as the 
; input image.

; ima_in: input image

; err_in: the error of flux in each pixel, should have the same shape as
; 	the input imag 

; dx: the x interval of the output image. The unit is 'pix'

; dy: the y interval of the output image. The unit is 'pix', could be 
; different with dx.

; ima_out: output image

; err_out: output error
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;code start;;;;;;;;;;;

ima_in=double(ima_in)
err_in=double(err_in)
dx=float(dx)
dy=float(dy)

var_in=err_in^2.

;;;;;;;;;;shape of input image;;;;;;;;;;
shape=size(ima_in)
nx_in=shape[1]
ny_in=shape[2]

;;;;;;;;;;;boundary of input image;;;;;
bdx_in=dindgen(nx_in+1)
bdy_in=dindgen(ny_in+1)

;;;;;;;;;;;shape of output image;;;;;;;;;
nx_out=floor(nx_in/dx)
ny_out=floor(ny_in/dy)


ima_out=dblarr(nx_out,ny_out)
var_out=dblarr(nx_out,ny_out)

;;;;;;;;;;;;;boundary of output image;;;;
bdx_out=dindgen(nx_out+1)*dx
bdy_out=dindgen(ny_out+1)*dy

;;;;;;;;;;position of the output boundary;;;;
kx=floor(bdx_out) < (nx_in-1) > 0
ky=floor(bdy_out) < (ny_in-1) > 0

;;;;;;;;;;find the flux in each output pixel;;
for i=0, nx_out-1 do begin
	for j=0, ny_out-1 do begin

		;;;;;x dir;;;;;;;
		ax=bdx_out[i]-bdx_in[kx[i]]
		bx=bdx_in[kx[i+1]+1]-bdx_out[i+1]

		;;;;;y dir;;;;;;;
		ay=bdy_out[j]-bdy_in[ky[j]]
		by=bdy_in[ky[j+1]+1]-bdy_out[j+1]

		;;;;calulate the flux in each pixel;;;;
		ima_out[i,j]=total(ima_in[kx[i]:kx[i+1],ky[j]:ky[j+1]]) - $
			(ax*total(ima_in[kx[i],ky[j]:ky[j+1]]) + bx*total(ima_in[kx[i+1],ky[j]:ky[j+1]]) + $
				ay*total(ima_in[kx[i]:kx[i+1],ky[j]]) + by*total(ima_in[kx[i]:kx[i+1],ky[j+1]])) + $
				(ax*ay*ima_in[kx[i],ky[j]] + ax*by*ima_in[kx[i],ky[j+1]] + $
					bx*ay*ima_in[kx[i+1],ky[j]] + bx*by*ima_in[kx[i+1],ky[j+1]])


		;;;;;calutate the variance in each piexl
		var_out[i,j]=total(var_in[kx[i]:kx[i+1],ky[j]:ky[j+1]]) - $
			(ax*total(var_in[kx[i],ky[j]:ky[j+1]]) + bx*total(var_in[kx[i+1],ky[j]:ky[j+1]]) + $
				ay*total(var_in[kx[i]:kx[i+1],ky[j]]) + by*total(var_in[kx[i]:kx[i+1],ky[j+1]])) + $
				(ax*ay*var_in[kx[i],ky[j]] + ax*by*var_in[kx[i],ky[j+1]] + $
					bx*ay*var_in[kx[i+1],ky[j]] + bx*by*var_in[kx[i+1],ky[j+1]])
	endfor
endfor

err_out=sqrt(var_out)

end




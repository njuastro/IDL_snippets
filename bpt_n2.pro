pro bpt_n2

;bpt diagram and different kind of maps(flux, rotation,gas velocity dispersion etc.) for IFU data.

;data set only for MaNAG MPL-6, you just need change the extension for other data release.

info = mrdfits('/Users/xiaoling/Desktop/project2/data/pro3/changing-look-agn.fits',1)
info1 = mrdfits('/Users/xiaoling/work/gal_samp_dap_v2_3_1.fits',1)
id_info1 = dblarr(n_elements(info))

n = n_elements(info)
for i = 0 , 0 do begin

    print,i

    plt = strtrim(string(info[i].plate),2)
    ifu = strtrim(string(info[i].ifu),2)
    id_info1[i] = where(plt eq info1.plate and ifu eq info1.ifudesign)
    ;stop
    dapdir = '/Users/xiaoling/work/MPL6-MAPS/'+plt+'/'+ifu+'/'
    dapname = 'manga-'+plt+'-'+ifu+'-MAPS-SPX-GAU-MILESHC.fits.gz'
    result_test = file_search(dapdir + dapname,count = num)
    if num eq 0 then continue
    snr = mrdfits(dapdir + dapname,5)
    fgas = mrdfits(dapdir + dapname,30)
    fgasivar = mrdfits(dapdir + dapname,31)
    fgasmask = mrdfits(dapdir + dapname,32)
    ew = mrdfits(dapdir + dapname,33)
    ew_ivar = mrdfits(dapdir + dapname,34)
    ew_mask = mrdfits(dapdir + dapname,35)
    v_gas = mrdfits(dapdir + dapname, 36)
    v_gas_ivar = mrdfits(dapdir + dapname, 37)
    v_gas_mask = mrdfits(dapdir + dapname, 38)
    siggas = mrdfits(dapdir + dapname,39)
    siggasivar = mrdfits(dapdir + dapname,40)
    siggasmask = mrdfits(dapdir + dapname,41)
    siginst = mrdfits(dapdir + dapname, 42)
    stellar_sigma = mrdfits(dapdir + dapname,18)
    stellar_sigma_ivar = mrdfits(dapdir + dapname,19)
    stellar_sigma_mask = mrdfits(dapdir + dapname,20)
    stellar_sigma_corr = mrdfits(dapdir + dapname,21)
    specindex = mrdfits(dapdir + dapname,44)
    specindex_ivar = mrdfits(dapdir + dapname,45)
    specindex_mask = mrdfits(dapdir + dapname,46)
    dn4000 = specindex[*,*,44]
    dn4000_ivar = specindex_ivar[*,*,44]
    dn4000_mask = specindex_mask[*,*,44]
    vstar = mrdfits(dapdir + dapname, 15)
    vstar_ivar = mrdfits(dapdir + dapname, 16)
    vstar_mask = mrdfits(dapdir + dapname, 17)
    ewha = ew[*,*,18]
    ewhaivar = ew_ivar[*,*,18]
    ewhamask = ew_mask[*,*,18]
    vha = v_gas[*,*,18]
    vo3 = v_gas[*,*,13]
    inst_ha = siginst[*,*,18]
    inst_o3 = siginst[*,*,13]
    sig_ha = siggas[*,*,18]
    sig_o3 = siggas[*,*,13]
    sigha = sqrt(sig_ha^2 - inst_ha^2)
    sigo3 = sqrt(sig_o3^2 - inst_o3^2)
    sigstar = sqrt(stellar_sigma^2 - stellar_sigma_corr^2)
    fhamask = fgasmask[*,*,18]
    fhaivar = fgasivar[*,*,18]
    fgaserr = 0 * fgasivar
    iok = where(fgasivar gt 0)
    if (iok[0] NE -1) THEN fgaserr[iok] = 1/(sqrt(fgasivar[iok]))
    
    o3 = fgas[*,*,13]
    o3err = fgaserr[*,*,13]
    hb = fgas[*,*,11]
    hberr = fgaserr[*,*,11]
    n2 = fgas[*,*,19]
    n2err = fgaserr[*,*,19]
    ha = fgas[*,*,18]
    haerr = fgaserr[*,*,18]
    

    s = size(ha)            
    nx = s[1]
    ny = s[2]
    xpos = dblarr(nx,ny)
    ypos= xpos
    xpos1 = dblarr(nx,ny)
    ypos1 = xpos1
    for xx = 0, nx-1 do begin 
       for yy = 0, ny-1 do begin
  
     xpos[xx,yy] = (xx-nx/2.)/2.
     ypos[xx,yy] = (yy-ny/2.)/2.
     xpos1[xx,yy] = xx
     ypos1[xx,yy] = yy
    endfor
      endfor
  
    bptax1 = alog10(n2/ha)
    bptay1 = alog10(o3/hb)
    indagn_center = where(((bptax1[nx/2,ny/2] lt 0.47) and (bptay1[nx/2,ny/2] ge (0.61/(bptax1[nx/2,ny/2]-0.47)+1.19))) or (bptax1[nx/2,ny/2] ge 0.47) and (o3[nx/2,ny/2]/o3err[nx/2,ny/2] ge 3. and hb[nx/2,ny/2]/hberr[nx/2,ny/2] ge 3 and n2[nx/2,ny/2]/n2err[nx/2,ny/2] ge 3 and ha[nx/2,ny/2]/haerr[nx/2,ny/2] ge 3) and abs(ewha[nx/2,ny/2]) gt 3)
    IF indagn_center eq -1 then continue
    indpix = where(fhamask eq 0 and fhaivar ne 0 and vstar_ivar ne 0 and snr ge 3)
    ;indpix = where(o3/o3err ge 3. and hb/hberr ge 3 and n2/n2err ge 3 and ha/haerr ge 3)
    if indpix[0] eq -1 then continue
 
vstar1 = 0 * vstar
vstar1[indpix] = vstar[indpix]
indpix_agn = where(((bptax1[indpix] lt 0.47) and (bptay1[indpix] ge (0.61/(bptax1[indpix]-0.47)+1.19))) or (bptax1[indpix] ge 0.47))
vstar1[indpix[indpix_agn]] = vo3[indpix[indpix_agn]]

;dfpsplot, '/Users/xiaoling/Desktop/project2/figure/pro3/changing-look-agn/'+strtrim(string(i+1),2)+'-'+plt+'-'+ifu+'.ps',/color
dfpsplot, '/Users/xiaoling/Desktop/project2/figure/pro3/changing-look-agn/'+strtrim(string(i+1),2)+'-'+plt+'-'+ifu+'.ps',/color
    device, decomposed=0
    sauron_colormap
    !P.MULTI = [0, 2, 5]
    !p.font=1
    !p.charsize=1.5
    !x.thick=1.2
    !y.thick=1.2

    imagename = '/Users/xiaoling/work/mpl6_rgb_image/'+plt+'-'+ifu+'.png'
    result = FILE_TEST(imagename)
    if result eq 1 then begin
    read_png, imagename, im
    imout = im[0:2,*,*]
    cgimage,imout, /keep
    endif

    ;dn4000 maps
    display_pixels, xpos[indpix], ypos[indpix], dn4000[indpix], range = [min(dn4000[indpix]),max(dn4000[indpix])],title = 'Dn 4000'+string(min(dn4000[indpix]))+string(max(dn4000[indpix])),pixelsize = 1,xtitle='  ',ytitle='  '
    cgCOLORBAR,NColors=255,minrange = min(dn4000[indpix]) , maxrange = max(dn4000[indpix]), /vertical, /right,font = 1,/normal,/fit,position = [0.88, 0.10, 0.90, 0.90]


    if min(ha[indpix]) ge 0 then minvalue = min(ha[indpix])
    if min(ha[indpix]) lt 0 then minvalue = 0

    ;ha flux maps
    display_pixels, xpos[indpix], ypos[indpix], ha[indpix], range = [minvalue,max(ha[indpix])],title = 'flux Halpha'+string(minvalue)+string(max(ha[indpix])),pixelsize = 1,xtitle='  ',ytitle=' '
    cgCOLORBAR,NColors=255,minrange = minvalue , maxrange = max(ha[indpix]), /vertical, /right,font = 1,/normal,/fit,position = [0.88, 0.10, 0.90, 0.90]
    cgcontour, ha, xpos, ypos,levels = [0.1, 0.5, 1.0,2.0, 4.0],c_colors = "white",/noerase,/overplot;,levels = [5, 10, 20, 30, 40, 50, 100, 200, 600],

    if min(o3[indpix]) ge 0 then minvalue = min(o3[indpix])
    if min(o3[indpix]) lt 0 then minvalue = 0

    ;oiii flux maps
    display_pixels, xpos[indpix], ypos[indpix], o3[indpix], range = [minvalue,max(o3[indpix])],title = 'flux OIII'+string(minvalue)+string(max(o3[indpix])),pixelsize = 1,xtitle='  ',ytitle=' '
    cgCOLORBAR,NColors=255,minrange = minvalue , maxrange = maxvalue, /vertical, /right,font = 1,/normal,/fit,position = [0.88, 0.10, 0.90, 0.90]

    if min(vstar[indpix]) gt -300 then minvalue = min(vstar[indpix])
    if min(vstar[indpix]) le -300 then minvalue = -300
    
    if max(vstar[indpix]) lt 300 then maxvalue = max(vstar[indpix])
    if max(vstar[indpix]) ge 300 then maxvalue = 300
    ; minvalue = -250
    ; maxvalue = 250
    display_pixels, xpos[indpix], ypos[indpix], vstar[indpix], title = 'V stars (DAP)'+string(minvalue)+string(maxvalue),pixelsize = 1,xtitle='  ',ytitle='  ',range = [minvalue,maxvalue]
    cgCOLORBAR,NColors=255,minrange = minvalue , maxrange = maxvalue, /vertical, /right,font = 1,/normal,/fit,position = [0.88, 0.10, 0.90, 0.90]


    if min(vha[indpix]) gt -250 then minvalue = min(vha[indpix])
    if min(vha[indpix]) le -250 then minvalue = -250
    
    if max(vha[indpix]) lt 250 then maxvalue = max(vha[indpix])
    if max(vha[indpix]) ge 250 then maxvalue = 250

    ; minvalue = -250
    ; maxvalue = 250

    ;ha rotation maps 
    display_pixels, xpos[indpix], ypos[indpix], vha[indpix], title = 'V Ha'+string(minvalue)+string(maxvalue),pixelsize = 1,xtitle='  ',ytitle='  ',range = [minvalue,maxvalue]
    cgCOLORBAR,NColors=255,minrange = minvalue , maxrange = maxvalue, /vertical, /right,font = 1,/normal,/fit,position = [0.88, 0.10, 0.90, 0.90]

    if min(sigstar[indpix]) gt 0 then minvalue = min(sigstar[indpix])
    if min(sigstar[indpix]) le 0 then minvalue = 0
    
    if max(sigstar[indpix]) lt 250 then maxvalue = max(sigstar[indpix])
    if max(sigstar[indpix]) ge 250 then maxvalue = 250

    ;stellar velocity dispersion maps
    display_pixels, xpos[indpix], ypos[indpix], sigstar[indpix], title = 'sigma stars (DAP)'+string(minvalue)+string(maxvalue),pixelsize = 1,xtitle='  ',ytitle='  ',range = [minvalue,maxvalue]
    cgCOLORBAR,NColors=255,minrange = minvalue , maxrange = maxvalue, /vertical, /right,font = 1,/normal,/fit,position = [0.88, 0.10, 0.90, 0.90]

    if min(sigha[indpix]) gt 0 then minvalue = min(sigha[indpix])
    if min(sigha[indpix]) le 0 then minvalue = 0
    
    if max(sigha[indpix]) lt 250 then maxvalue = max(sigha[indpix])
    if max(sigha[indpix]) ge 250 then maxvalue = 250
    
    ;gas velocity dispersion maps
    display_pixels, xpos[indpix], ypos[indpix], sigha[indpix], title = 'sigma gas'+string(minvalue)+string(maxvalue),pixelsize = 1,xtitle='  ',ytitle='  ',range = [minvalue,maxvalue]
    cgCOLORBAR,NColors=255,minrange = minvalue , maxrange = maxvalue, /vertical, /right,font = 1,/normal,/fit,position = [0.88, 0.10, 0.90, 0.90]

    ;indplot = where(o3/o3err ge 3. and hb/hberr ge 3 and n2/n2err ge 3 and ha/haerr ge 3)
    indplot = where(fhamask eq 0 and fhaivar ne 0 and vstar_ivar ne 0 and snr ge 3)
    if indplot[0] ne -1 then begin
    npix = n_elements(indplot)
    flag = replicate(100.,npix)

    ;separate SF galaxies and AGNs using BPT (see Kauffmann et al. 2003)
    bptax = alog10(n2[indplot]/ha[indplot])
    bptay = alog10(o3[indplot]/hb[indplot])
    indsf = where((bptax lt 0.05) and (bptay le (0.61/(bptax-0.05)+1.3)))
    indcomp = where((bptay ge (0.61/(bptax-0.05)+1.3)) and (bptay le (0.61/(bptax-0.47)+1.19)))
    indagn = where(((bptax lt 0.47) and (bptay ge (0.61/(bptax-0.47)+1.19))) or (bptax ge 0.47),ct)
    indliner = where((((bptax lt 0.47) and (bptay ge (0.61/(bptax-0.47)+1.19))) or (bptax ge 0.47)) and (bptax ge alog10(0.6) and bptay le alog10(3.)),ct)
    indsey = where(bptax ge alog10(0.6) and bptay ge alog10(3.))
    
    plot,bptax,bptay,psym=cgsymcat(16),symsize = 0.5,xr=[-1.5,0.5],yr=[-1.2,1.5],/xstyle,/ystyle,title='BPT diagram',$
    ytitle=textoidl('log [OIII]\lambda5007/H\beta'),xtitle=textoidl('log [NII]\lambda6583/H\alpha'),$
    charsize = 1.5
    
    if indagn[0] ne -1 then begin
    flag[indagn] = 1.5 
    oplot,bptax[indagn],bptay[indagn],color=djs_icolor('red'),psym=cgsymcat(16),symsize = 0.5    
    endif

    if indsf[0] ne -1 then begin
    flag[indsf] = - 1.5 
    oplot,bptax[indsf],bptay[indsf],color=djs_icolor('blue'),psym=cgsymcat(16),symsize = 0.5    
    endif
    if indsey[0] ne -1 then begin    
    flag[indsey] = 1.5
    oplot,bptax[indsey],bptay[indsey],color=djs_icolor('red'),psym=cgsymcat(16),symsize = 0.5
    endif
    if indliner[0] ne -1 then begin
    flag[indliner] = 0.5
    oplot,bptax[indliner],bptay[indliner],color=djs_icolor('yellow'),psym=cgsymcat(16),symsize = 0.5    
    endif
    if indcomp[0] ne -1 then begin
    flag[indcomp] = 0.
    oplot,bptax[indcomp],bptay[indcomp],color=djs_icolor('green'),psym=cgsymcat(16),symsize = 0.5    
    endif
    xxx = (indgen(20)-15.)*0.1
    xxx1 = (indgen(20)-19.5)*0.1
    yyy1 = 0.61/(xxx1-0.05)+1.3
    yyy2 = 0.61/(xxx-0.47)+1.19
    oplot,xxx1,yyy1
    oplot,xxx,yyy2,linestyle=2
    oplot,[alog10(0.6),alog10(0.6)],[-100,100]
    oplot,[alog10(0.6),100],[alog10(3.),alog10(3.)]
    display_pixels, xpos[indplot], ypos[indplot], flag*120,$ 
           range = [-250,250], title = strtrim(string(plt),2)+'-'$
             +strtrim(ifu,2)+'Ka03',pixelsize = 1,xtitle='  ',ytitle='  ',charsize = 1.5,charthick=2
    endif
dfpsclose
endfor
;MWRFITS,info1[id_info1],'/Users/xiaoling/Desktop/project2/data/pro3/changing-look-agn-info.fits'
stop
end

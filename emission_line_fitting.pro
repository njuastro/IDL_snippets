PRO emission_line_fitting

       ;purpose:fitting ha and hb region and obtain sigma，flux，error etc.


            



            
       
            stack_fits = file_search('/Users/xiaoling/Desktop/work/ppxf_stack_v4/stack_ppxf_fits','*ppxf.fits',count = num)
        


            fit_plate = strarr(num)
            fit_ifu = fit_plate
            fit_sigma = dblarr(num)
            fit_sigma_err = fit_sigma
            fit_flux = fit_sigma
            fit_flux_err = fit_sigma

            for i = 0,num - 1 do begin
            split = strsplit(stack_fits[i],'-',/extract)
            plt = split[1]
            ifu = split[2]
            spectral = mrdfits(stack_fits[i],1)
            ind_wave = where(spectral.wave ge 6530. and spectral.wave le 6610.)
            wave = spectral[ind_wave].wave
            flux = spectral[ind_wave].fit_diff
            flux_err = spectral[ind_wave].flux_err
            indwave_NII_1 = where(spectral.wave ge 6538. and spectral.wave le 6558.)
            flux_NII_1 = total(spectral[indwave_NII_1].fit_diff)
            indwave_Ha = where(spectral.wave ge 6552. and spectral.wave le 6572.)
            flux_Ha = total(spectral[indwave_Ha].fit_diff)
            indwave_NII_2 = where(spectral.wave ge 6573. and spectral.wave le 6593.)
            flux_NII_2 = total(spectral[indwave_NII_2].fit_diff)



            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;           
            ;;;;fit Halpha region, including Halpha;;;;
            ;;;;NII [6549,6583];;;;;;;;;;;;;;;;;;;;;
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            pi = replicate({limited:[0,0],limits:[0.,0.],mpside:2},9)
            
            pi[3].limited=[1,1]
            pi[3].limits=[6552.,6572.]
            pi[4].limited=[1,1]
            pi[4].limits=[0.5D,300.D]
            pi[5].limited=[1,1]
            pi[5].limits=[0.1D,flux_Ha+30000]

            expr = 'gauss1(x,p[0:2]) + gauss1(x,p[3:5]) + gauss1(x,p[6:8])'
            start = [6548., 2., flux_NII_1, 6562., 2., flux_Ha, 6583., 2, flux_NII_2]
            poreg2 = mpfitexpr(expr, wave, flux, flux_err, start, parinfo = pi, perror = pe, bestnorm=ch,/quiet)
            
            ind_wave1 = where(spectral.wave ge 4800. and spectral.wave le 5100.)
            wave1 = spectral[ind_wave1].wave
            flux1 = spectral[ind_wave1].fit_diff
            flux_err1 = spectral[ind_wave1].flux_err
            indwave_Hb = where(spectral.wave ge 4851. and spectral.wave le 4871.)
            flux_Hb = total(spectral[indwave_Hb].fit_diff)
            indwave_OIII_1 = where(spectral.wave ge 4948. and spectral.wave le 4968.)
            flux_OIII_1 = total(spectral[indwave_OIII_1].fit_diff)
            indwave_OIII_2 = where(spectral.wave ge 4997. and spectral.wave le 5017.)
            flux_OIII_2 = total(spectral[indwave_OIII_2].fit_diff)
            pi1 = replicate({limited:[0,0],limits:[0.,0.],mpside:2},9)




            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;           
            ;;;;fit Hbeta region, including Hbeta;;;;
            ;;;;OIII [4959,5007];;;;;;;;;;;;;;;;;;;;;
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            
            pi1[0].limited=[1,1]
            pi1[0].limits=[4851.,4871.]
            pi1[1].limited=[1,1]
            pi1[1].limits=[0.5D,300.D]
            pi1[2].limited=[1,1]
            pi1[2].limits=[0.1D,flux_Hb+30000]

            expr1 = 'gauss1(x,p[0:2]) + gauss1(x,p[3:5]) + gauss1(x,p[6:8])'
            start1 = [4861., 2., flux_Hb, 4959., 2., flux_OIII_1, 5007., 2, flux_OIII_2]
            poreg3 = mpfitexpr(expr1, wave1, flux1, flux_err1, start1, parinfo = pi1, perror = pe1, bestnorm=ch, /WEIGHTS,/quiet)




            fit_plate[i] = plt
            fit_ifu[i] = ifu
            fit_sigma[i] = ((poreg2[4])/(6562.8)) *(3 *10.^5)
            fit_sigma_err[i] = (pe[4]/(6562.8)) *(3 *10.^5)
            fit_flux[i] = poreg2[5]
            fit_flux_err[i] = pe[5]

;stop
xxss=12
yyss=12
set_plot,'ps'
LOADCT,39
device, filename='/Users/xiaoling/Desktop/gfit_figure/'+ 'manga'+'-'+plt+'-'+ifu+'-'+'gauss.ps',xsize=35,ysize=10,/ENCAPSULATED,/color

!p.font=0
!p.charsize=2.5
!x.thick=2.0
!y.thick=2.0
!p.multi = [0,3,1,0,0]

cgplot,[0],[0],/xs,/ys,/nodata,xr = [6520,6620],yr = [-0.2,max(flux)+0.2],title = 'manga'+'-'+plt+'-'+ifu,xtitle = textoidl('wavelength ($\Angstrom$)'),ytitle =textoidl('Flux (10^{-17} erg s^{-1} cm^{-2} $\Angstrom$^{-1})')         
cgplot, wave, flux,thick = 5.0,/overplot 
cgplot,wave,gauss1(wave,poreg2[0:2]) + gauss1(wave,poreg2[3:5]) + gauss1(wave,poreg2[6:8]),color=djs_icolor('red'),thick = 5.0,/overplot
;cgplot,[6562.8,6562.8],[-20000,20000],color=djs_icolor('green'),thick = 5.0,/overplot


cgplot,[0],[0],/xs,/ys,/nodata,xr = [4800,5100],yr = [-20,max(flux1)+20],title = 'manga'+'-'+plt+'-'+ifu,xtitle = textoidl('wavelength ($\Angstrom$)'),ytitle =textoidl('Flux (10^{-17} erg s^{-1} cm^{-2} $\Angstrom$^{-1})')         
cgplot, wave1, flux1,thick = 5.0,/overplot
cgplot,wave1,gauss1(wave1,poreg3[0:2]) + gauss1(wave1,poreg3[3:5]) + gauss1(wave1,poreg3[6:8]),color=djs_icolor('red'),thick = 5.0,/overplot

device,/close          

out_put_result = replicate({WAVE:double(0.),FLUX:double(0.),MODEL_FLUX:double(0.)},n_elements(wave))
out_put_result.WAVE = wave
out_put_result.FLUX = flux
out_put_result.MODEL_FLUX = gauss1(wave,poreg2[0:2]) + gauss1(wave,poreg2[3:5]) + gauss1(wave,poreg2[6:8])
MWRFITS,out_put_result,'/Users/xiaoling/Desktop/ha/'+ 'manga'+'-'+plt+'-'+ifu+'-'+'gauss.fits'



            ENDFOR

           fit_gauss_result = replicate({PLATE:'',IFU:'',SIGMA:double(0.),ERR_SIGMA:double(0.),FLUX:double(0.),ERR_FLUX:double(0.)},n_elements(stack_fits))
           fit_gauss_result.PLATE = fit_plate
           fit_gauss_result.IFU = fit_ifu
           fit_gauss_result.SIGMA = fit_sigma
           fit_gauss_result.ERR_SIGMA = fit_sigma_err
           fit_gauss_result.FLUX = fit_flux
           fit_gauss_result.ERR_FLUX = fit_flux_err


           ;MWRFITS,fit_gauss_result,'/Users/xiaoling/Desktop/work/test_corr_uncorr/gauss_fit_data_un/'+ 'mangafit'+'_gauss_ha.fits'





STOP
END
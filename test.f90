program test
    use fp_sim
    use rsplot
    real :: reflect_1 = 0.95, &
            reflect_2 = 0.98
    double precision :: length=0.3, increment=5e-9, x, y
    double precision, dimension(:,:), allocatable :: spectrum_in, spectrum_out
    double precision :: temp=6000., noise=0.005, min=600e-9, max=1100e-9, resolution=0.1e-9
    double precision :: R1=0.5, R2=0.5, length_init = 0.01, start_lamda = 600e-9, finish_lamda = 1050e-9, step = 1e-9
    double precision, dimension(6,3) :: peaks
    integer :: counter = 0, counter2
    real :: upper_in, lower_in, lower_out, upper_out
    real, dimension(:), allocatable :: spect_x_in, spect_y_in, spect_x_out, spect_y_out
    double precision :: end_point, start_point

!===========================================================================================
!    Create spectrum for input
!      
!===========================================================================================
     
     peaks(1,:) = (/627d-9, 2d-9, 1d-1/)
     peaks(2,:) = (/656d-9, 2d-9,  2d-1/)
     peaks(3,:) = (/687d-9, 3d-9, 3d-1/)
     peaks(4,:) = (/759d-9, 2d-9, 2d-1/)
     peaks(5,:) = (/823d-9, 3d-9, 1d-1/)
     peaks(6,:) = (/898d-9, 4d-9, 1d-1/) 
     counter = nint((max-min)/resolution)
     allocate (spectrum_in(2,counter))
     allocate (spect_x_in(counter))
     allocate (spect_y_in(counter))
     spectrum_in = create_spectrum(temp=temp, peaks=peaks, noise=noise, min=min, max=max, resolution=resolution)
     spect_x_in = real(spectrum_in(1,:))
     spect_y_in = real(spectrum_in(2,:))
     upper_in = maxval(spect_y_in)
     lower_in = minval(spect_y_in)
!===========================================================================================
!    Simulate Fabry-Perot Spectrometer
!
!===========================================================================================
     start_point = nint((length_init/start_lamda))*start_lamda
     end_point = start_point + (finish_lamda - start_lamda)                                                
     counter2 = nint((end_point - start_point)/step)
     allocate (spectrum_out(3,counter2))
     allocate (spect_x_out(counter2))
     allocate (spect_y_out(counter2))
     spectrum_out = scan_cavity(spectrum_in, R1, R2, length_init, start_lamda, finish_lamda, step)
     spect_x_out = real(spectrum_out(2,:))
     spect_y_out = real(spectrum_out(3,:))
     upper_out = maxval(spect_y_out)
     lower_out = minval(spect_y_out)
!============================================================================================
!    Plot the two spectra
!
!============================================================================================    
     call ps_init(pform='A0_P', rot90=0, fileout='spectrum')
     call ps_frame(1, real(min), real(lower_in), real(max),upper_in*1.1, &
                   xbl=40.,ybl=650. , clip=1, &
                   fill=ps_color([.94,.94,.94]))
     call ps_axis(1, ax='Xx', glp=0, d_tm=5., title= &
          'Wavelength [m]')
     call ps_axis(1, ax='Yy', title='Intensity', glp=0)
     call ps_plot(1, counter, spect_x_in, spect_y_in, ci=1)
     call ps_frame(2, real(start_lamda), real(lower_out), real(finish_lamda), &
                   upper_out, xbl=40., &
                   ybl = 40., clip=1, fill=ps_color([.94,.94,.94]))
     call ps_axis(2, ax='Xx', glp=0, d_tm=5., title= 'Wavelength [m]')     
     call ps_axis(2, ax='Yy', title='Intensity', glp=0) 
     call ps_plot(2, counter2, spect_x_out, spect_y_out, ci=1)    
     call ps_exit    
     
end program test

$(document).ready(function () {
    var var_sidebar = `<nav id="sidebar"> 
         
         
                    <div class="sidebar-header"> 
                        <a href="/website/html/" class="sidelink"> 
                            <img class="img-fluid" src="/website/images/logo.png" alt="Cornell Logo"> 
                        </a> 
                    </div> 
         
                    <ul class="list-unstyled components"> 
                        <li id="sidehome"> 
                            <a href="/website/html/" class="sidelink">Home</a> 
                        </li> 
                        <li id="sidebg"> 
                            <a href="#bgSubmenu" data-toggle="collapse" aria-expanded="false" class="sidelink">Background</a> 
                            <ul class="collapse list-unstyled" id="bgSubmenu" data-parent="#sidebar"> 
                                <li id="sideanalysis"> 
                                    <a href="/website/html/background/analysis_on_fractals.html" class="sidelink">Analysis on Fractals</a> 
                                </li> 
                                <li id="sidepoly"> 
                                    <a href="/website/html/background/polynomials_on_sg.html" class="sidelink">Polynomials on SG</a> 
                                </li> 
                            </ul> 
                        </li> 
         
                        <li id="sidetheory"> 
                            <a href="#theorySubmenu" data-toggle="collapse" aria-expanded="false" class="sidelink">Theory</a> 
                            <ul class="collapse list-unstyled" id="theorySubmenu" data-parent="#sidebar"> 
                                <li id="sideSOP"> 
                                    <a href="/website/html/theory/sobolev_ops.html" class="sidelink">Sobolev Orthogonal Polynomials</a> 
                                </li> 
                                <li id="sideCP"> 
                                    <a href="/website/html/theory/chebyshev_ops.html" class="sidelink">Chebyshev Polynomials</a> 
                                </li> 
                            </ul> 
                        </li> 
                        <li id="sidenumerical"> 
                            <a href="#numSubmenu" data-toggle="collapse" aria-expanded="false" class="sidelink">Numerical Results</a> 
                            <ul class="collapse list-unstyled" id="numSubmenu" data-parent="#sidebar"> 
                                <li id="sideplot"> 
                                    <a href="/website/html/numerical_results/plots.html" class="sidelink">Plots</a> 
                                </li> 
                                <li id="sidecomp"> 
                                    <a href="/website/html/numerical_results/comparisons.html" class="sidelink">Comparisons</a> 
                                </li> 
                                <li id="sidezero"> 
                                    <a href="/website/html/numerical_results/zeros.html" class="sidelink">Zeros</a> 
                                </li> 
                                <li id="sideapp"> 
                                        <a href="/website/html/numerical_results/applications.html" class="sidelink">Applications</a> 
                                </li> 
                            </ul> 
                        </li> 
                        <li id="siderec"> 
                            <a href="#resourceSubmenu" data-toggle="collapse" aria-expanded="false" class="sidelink">Resources</a> 
                            <ul class="collapse list-unstyled" id="resourceSubmenu" data-parent="#sidebar"> 
                                <li id="sidecode"> 
                                    <a href="/website/html/resources/programs.html" class="sidelink">Software and Programs</a> 
                                </li> 
                                <li id="sidepresent"> 
                                    <a href="/website/html/resources/presentations.html" class="sidelink">Presentations</a> 
                                </li> 
                                <li id="sideref"> 
                                    <a href="/website/html/resources/references.html" class="sidelink">References</a> 
                                </li> 
                            </ul> 
                        </li> 
                        <li id="sidecontact"> 
                            <a href="/website/html/contact.html">Contact Us</a> 
                        </li> 
                    </ul> 
         
                     
                </nav>`;


        var var_imgModal = `<div id="myModal" class="modal"> 
         
        <span class="close"><i class="fas fa-window-close" aria-hidden="true"></i></span> 
 
        <img class="modal-content" id="img01"> 
 
        <div id="caption"></div> 
        </div>`;

        var var_scroll = `<a id="back-to-top" href="#" class="btn btn-dark btn-lg back-to-top" role="button"><i class="fas fa-chevron-up"></i></a>`;
        var var_overlay = `<div class="overlay"></div>`

        var var_button = `<button type="button" id="sidebarCollapse" class="btn btn-info pull-left"  style="margin-top:5px;">
                                Menu <i class="fas fa-align-justify"></i>
                            </button>`;

    $(var_sidebar).prependTo($("#main"));
    $(var_imgModal).prependTo($("#main"));
    $(var_scroll).prependTo($("#main"));
    $(var_overlay).prependTo($("#main"));
    $(var_button).prependTo($("#content"));

    $("#sidebar").mCustomScrollbar({
        theme: "minimal"

    });

    $('#dismiss, .overlay').on('click', function () {
        $('#sidebar').removeClass('active');
        $('.overlay').removeClass('active');
    });

    $('#sidebarCollapse').on('click', function () {
        $('#sidebar').addClass('active');
        $('.overlay').addClass('active');
        $('.collapse.in').toggleClass('in');
        $('a[aria-expanded=true]').attr('aria-expanded', 'false');
    });

    $('.imgModal').on('click', function () {
        var modal = document.getElementById("myModal");
        var modalImg = document.getElementById("img01");
        var captionText = document.getElementById("caption");
        modal.style.display = "block";
        modalImg.src = this.src;
        captionText.innerHTML = this.alt;
    });

    $('.close').on('click', function () {
        var modal = document.getElementById("myModal");
        modal.style.display = "none";
    })

    $(window).on('click', function (e) {
        var modal = document.getElementById("myModal");
        if (event.target == modal) {
            modal.style.display = "none";
        }
    })

    $(window).scroll(function () {
        if ($(this).scrollTop() > 50) {
            $('#back-to-top').fadeIn();
        } else {
            $('#back-to-top').fadeOut();
        }
    });
    // scroll body to 0px on click
    $('#back-to-top').click(function () {

        $('body,html').animate({
            scrollTop: 0
        }, 400);
        return false;
    });
});   

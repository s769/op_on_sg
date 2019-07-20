$(document).ready(function () {
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

    $('.imgModal').on('click', function(){
        console.log('sus');
        var modal = document.getElementById("myModal");
        var modalImg = document.getElementById("img01");
        var captionText = document.getElementById("caption");
        modal.style.display = "block";
        modalImg.src = this.src;
        captionText.innerHTML = this.alt;
    });

    $('.close').on('click', function(){
        var modal = document.getElementById("myModal");
        modal.style.display = "none";
    })
    
    $(window).on('click', function(e){
        var modal = document.getElementById("myModal");
        if (event.target == modal){
            modal.style.display = "none";
        }
    })

});   
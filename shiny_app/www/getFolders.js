document.getElementById("fileIn").addEventListener("change", function(e) {

        let files = e.target.files;
        var relativePath = files[0].webkitRelativePath;
        var folder = relativePath.split("/");

        Shiny.onInputChange("mydata", folder[0]);

});
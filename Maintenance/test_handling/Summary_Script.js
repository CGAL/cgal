const $mainContainer = $('#main_container');
const $platformContainer = $('#platform_container');
const $packageContainer = $('#package_container');
let release = "";
let packages = [];

function getAllTestDirectories(data) {
    return data.platforms.flatMap(platform => platform.test_directories.map(directory => directory.test_directory));
}

function clearPackagesOptions() {
    $('#packagesSelector option').each(function() {
        if ($(this).val() !== 'all' && !packages.includes($(this).val())) {
            $(this).remove();
        }
    });
}

function filterByPlatform(platform) {
    const $packageSelector = $('#packagesSelector');
    $packageSelector.prop('disabled', platform !== 'all');
    $platformContainer.find('.platform').each(function() {
        const $this = $(this);
        if (platform === 'all' || $this.hasClass(platform)) {
            $this.show();
        } else {
            $this.hide();
        }
    });
}

function filterByPackage(package) {
    const $platformSelector = $('#platformSelector');
    const $packageContainers = $packageContainer.find('.package');
    if (package === 'all') {
        $platformSelector.prop('disabled', false);
        $platformContainer.show();
        $packageContainer.hide();
    } else {
        $platformSelector.prop('disabled', true);
        $platformContainer.hide();
        $packageContainer.show();
    }
    $packageContainers.each(function() {
        const $this = $(this);
        if (package === 'all' || $this.hasClass(package)) {
            $this.show();
        } else {
            $this.hide();
        }
    });
}

function filterByLetter(letter) {
    const $letterContainers = $mainContainer.find('.letter_container');
    $letterContainers.each(function() {
        const $this = $(this);
        if (letter === 'all' || $this.hasClass(letter)) {
            $this.show();
        } else {
            $this.hide();
        }
        if ($this.children().length <= 2) {
            $this.hide();
        }
    });
}

function search() {
    const searchTerm = $('#searchInput').val().toLowerCase();
    const releaseType = $('#releaseSelector').val();
    const $resultsContainer = $('#searchResults');
    $resultsContainer.empty();

    if (!searchTerm) {
        $resultsContainer.append('<p>Please enter a search term.</p>');
        return;
    }

    const matchingDirectories = [];

    window.data.platforms.forEach(platform => {
        platform.test_directories.forEach(directory => {
            if (directory.content.toLowerCase().includes(searchTerm)) {
                matchingDirectories.push({
                    test_directory: directory.test_directory,
                    platform_name: platform.name
                });
            }
        });
    });

    if (matchingDirectories.length === 0) {
        $resultsContainer.append('<p>No matching directories found.</p>');
    } else {
        matchingDirectories.forEach(match => {
            const link = `${window.data.release}/${match.test_directory}/TestReport_${match.platform_name}.gz`;
            $resultsContainer.append(`<p><a href="${link}" target="_blank">${match.platform_name} - ${match.test_directory} - <strong>${window.data.release}</strong></a></p>`);
        });
    }
}

function packageContainer(platforms) {
    const testDirectories = {};

    platforms.forEach(platform => {
        platform.test_directories.forEach(testDir => {
            if (!testDirectories[testDir.test_directory]) {
                testDirectories[testDir.test_directory] = {};
            }
            if (!testDirectories[testDir.test_directory][testDir.letters]) {
                testDirectories[testDir.test_directory][testDir.letters] = [];
            }
            testDirectories[testDir.test_directory][testDir.letters].push({
                platformName: platform.name,
                content: testDir.content
            });
        });
    });

    for (const [testDirectory, letters] of Object.entries(testDirectories)) {
        const $container = $('<div>', { class: 'package ' + testDirectory }).appendTo($packageContainer);
        $('<h2>').text(testDirectory).appendTo($container);

        for (const [letter, platformDetails] of Object.entries(letters)) {
            const $letterContainer = $('<div>', {
                class: 'letter_container ' + letter,
            }).appendTo($container);
            $('<h3>').text(letter).appendTo($letterContainer);

            platformDetails.forEach(detail => {
                const { platformName, content } = detail;

                const $platformContainer = $('<div>', { class: 'platform-container' }).appendTo($letterContainer);

                const $platformLink = $('<a>', {
                    href: `${release}/${testDirectory}/TestReport_${platformName}.gz`,
                    text: platformName,
                    class: 'platform-link'
                }).appendTo($platformContainer);

                const $contentSpan = $('<pre>', {
                    text: content,
                    class: 'summary-content',
                    css: { display: 'none' }
                }).appendTo($letterContainer);

                if (content.length > 0) {
                    const $toggleButton = $('<button>', {
                        class: 'toggle-button',
                        text: 'Show More',
                        click: function() {
                            if ($contentSpan.is(':hidden')) {
                                $contentSpan.show().css('background-color', '#D0D0E0');
                                $(this).text('Show Less');
                            } else {
                                $contentSpan.hide().css('background-color', 'transparent');
                                $(this).text('Show More');
                            }
                        }
                    }).appendTo($platformContainer);
                }

                $('<br>').appendTo($letterContainer);
            });
        }
    }
}

function platformContainer(platforms) {
    platforms.forEach(platform => {
        const $container = $('<div>', { class: 'platform ' + platform.name }).appendTo($platformContainer);
        $container.html(`<h2>Results of ${platform.name}</h2>`);
        const tplArray = platform.tpl;
        const $toggleButton = $('<button>', {
            text: 'Third Party Libraries',
            class: 'tpl-toggle-button toggle-button',
            click: function() {
                $tplTable.toggle();
            }
        }).appendTo($container);
        const $tplTable = $('<table>', {
            class: 'tablesorter',
            css: {
            display: 'none',
            maxWidth: '300px'
            }
        }).appendTo($container);
        const $thead = $('<thead>').appendTo($tplTable);
        const $tbody = $('<tbody>').appendTo($tplTable);
        const $headerRow = $('<tr>');
        $('<th>', { text: 'Library' }).appendTo($headerRow);
        $('<th>', { text: 'Version' }).appendTo($headerRow);
        $headerRow.appendTo($thead);
        tplArray.forEach(tpl => {
            const $row = $('<tr>').append(
                $('<td>').text(tpl.name),
                $('<td>').text(tpl.version || 'N/A'),
            ).appendTo($tbody);
            $row.addClass('tpl-row');
            $row.click(function() {
                showVersionsForTPL(tpl.name);
            });
        });
        const letters = ['n', 'w', 'o', 'r'];
        letters.forEach(letter => {
            const $letterContainer = $('<div>', { class: 'letter_container ' + letter }).appendTo($container);
            $('<h3>').text(letter).appendTo($letterContainer);
            const testDirectoriesForLetter = platform.test_directories.filter(directory => directory.letters === letter);
            testDirectoriesForLetter.forEach(directory => {
                const $directoryContainer = $('<div>', { class: 'directory_container' }).appendTo($letterContainer);
                const $directoryName = $('<a>', {
                    href: `${release}/${directory.test_directory}/TestReport_${platform.name}.gz`,
                    text: `${directory.test_directory}  `
                }).appendTo($directoryContainer);
                const $contentSpan = $('<pre>', {
                    text: directory.content,
                    class: 'summary-content',
                    css: { display: 'none' }
                }).appendTo($letterContainer);
                if (directory.content.length > 0) {
                    const $toggleButton = $('<button>', {
                        class: 'toggle-button',
                        text: 'Show More',
                        click: function() {
                            if ($contentSpan.is(':hidden')) {
                                $contentSpan.show().css('background-color', '#D0D0E0');
                                $(this).text('Show Less');
                            } else {
                                $contentSpan.hide().css('background-color', 'transparent');
                                $(this).text('Show More');
                            }
                        }
                    }).appendTo($directoryContainer);
                    $('<br>').appendTo($letterContainer);
                }
            });
            if ($letterContainer.children().length <= 2) {
                $letterContainer.hide();
            }
        });
    });
}

function openAll() {
    $('.summary-content').show().css('background-color', '#D0D0E0');
    $('.toggle-button').text('Show Less');
}

function closeAll() {
    $('.summary-content').hide().css('background-color', 'transparent');
    $('.toggle-button').text('Show More');
}

function showVersionsForTPL(tplName) {
    const $modal = $('#tplModal');
    const $modalTitle = $('#tplModalTitle');
    const $modalTable = $modal.find('table');
    const $modalBody = $modalTable.find('tbody');
    $modalBody.empty();
    $modalTitle.text(`Versions of ${tplName} across platforms`);
    let tplFound = false;
    window.data.platforms.forEach(platform => {
        const matchingTPL = platform.tpl.find(tpl => tpl.name === tplName);
        if (matchingTPL) {
            tplFound = true;
            $modalBody.append(`
                <tr class="modal-table-row">
                    <td>${platform.name}</td>
                    <td>${matchingTPL.version || 'N/A'}</td>
                </tr>
            `);
        }
    });
    if (!tplFound) {
        $modalBody.append('<tr><td colspan="3">No versions of this TPL found across platforms.</td></tr>');
    }
    $modalTable.trigger("update");
    $modal.show();
    $('.close').click(function() {
        $modal.hide();
    });
    $(window).click(function(event) {
        if (event.target == $modal[0]) {
            $modal.hide();
        }
    });
}

function main() {
    const url = searchURLs["current"];
    $.getJSON(url, data => {
        window.data = data;
        release = data.release;
        packages = getAllTestDirectories(data);
        clearPackagesOptions();
        platformContainer(data.platforms);
        packageContainer(data.platforms);
        $packageContainer.hide();
        $(document).ready(function() {
            $("table.tablesorter").tablesorter();
        });
        const urlParams = new URLSearchParams(window.location.search);
        const platform = urlParams.get('platform');
        if (platform) {
            console.log("Platform is " + platform);
            filterByPlatform(platform);
        }
        const package = urlParams.get('package');
        if (package) {
            console.log("Package is " + package);
            filterByPackage(package);
        }
        $('#open-all').click(openAll);
        $('#close-all').click(closeAll);
        $('#search-button').click(search);
    }).fail(error => console.log(error));
}

$(document).ready(main);
